#ifndef PACKET_HPP
#define PACKET_HPP

#include "header.hpp"
#include <stdlib.h>
#include <vector>

struct packet {
    // unsigned int type;
    unsigned int mode_number;
    unsigned int mapping_number;
    int n;
    bool previous_window_flag;
    bool next_window_flag;
    vector<double> window;
    int window_center;
    int left_window_start;
    int left_window_end;
    int left_n;
    int right_window_start;
    int right_window_end;
    int right_n;

    vector<vector<double>> floor_output;
    vector<bool> no_residue;
    vector<vector<double>> r_output;
    vector<vector<double>> r_output_tmp;

    vector<vector<double>> previous_y;
    vector<vector<double>> Y;
    vector<vector<double>> Y_ret;

    packet(identification &id) {
        n = max(id.blocksize_0, id.blocksize_1);
        previous_y.resize(id.audio_channels);
        Y.resize(id.audio_channels);
        Y_ret.resize(id.audio_channels);
        for (int i = 0; i < id.audio_channels; ++i)
            previous_y[i].assign(n, 0);
    }

    void decode_window(io_buf &in, identification &id, setup &s, modes &mode, mappings &mapping) {
        unsigned int blockflag = mode.blockflag;
        previous_window_flag = next_window_flag = 0;
        n = (blockflag == 0) ? (int)id.blocksize_0 : (int)id.blocksize_1;
        if (blockflag == 1) {
            previous_window_flag = in.read_u(1);
            next_window_flag = in.read_u(1);
        }

        window.clear();
        window.resize(n);
        window_center = n / 2;
        if (blockflag == 1 && previous_window_flag == 0) {
            left_window_start = n / 4 - id.blocksize_0 / 4;
            left_window_end = n / 4 + id.blocksize_0 / 4;
            left_n = id.blocksize_0 / 2;
        } else {
            left_window_start = 0;
            left_window_end = window_center;
            left_n = n / 2;
        }
        if (blockflag == 1 && next_window_flag == 0) {
            right_window_start = n * 3 / 4 - id.blocksize_0 / 4;
            right_window_end = n * 3 / 4 + id.blocksize_0 / 4;
            right_n = id.blocksize_0 / 2;
        } else {
            right_window_start = window_center;
            right_window_end = n;
            right_n = n / 2;
        }
        for (int i = 0; i < left_window_start; ++i)
            window[i] = 0;
        for (int i = left_window_start; i < left_window_end; ++i)
            window[i] = sin(M_PI / 2 * pow(sin((i - left_window_start + 0.5) / left_n * (M_PI / 2)), 2));
        for (int i = left_window_end; i < right_window_start; ++i)
            window[i] = 1;
        for (int i = right_window_start; i < right_window_end; ++i)
            window[i] = sin(M_PI / 2 * pow(sin((i - right_window_start + 0.5) / right_n * (M_PI / 2) + M_PI / 2), 2));
        for (int i = right_window_end; i < n; ++i)
            window[i] = 0;
    }

    void decode_floor_curve(io_buf &in, identification &id, setup &s, modes &mode, mappings &mapping) {
        // floor curve decode
        // vector<vector<double>> floor_output;
        floor_output.clear();
        floor_output.resize(id.audio_channels);
        for (unsigned int t = 0; t < id.audio_channels; ++t)
            floor_output[t].assign(n / 2, 0);

        no_residue.clear();
        no_residue.resize(id.audio_channels);
        for (int i = 0; i < id.audio_channels; ++i)
            no_residue[i] = 0;

        for (unsigned int t = 0; t < id.audio_channels; ++t) {
            /*
               unsigned int submap_number;
               if(mapping.submap > 1) submap_number = mapping.mux[t];
               else submap_number = 0;
             */
            unsigned int submap_number = mapping.mux[t];
            unsigned int floor_number = mapping.submap_floor[submap_number];
            auto &fl = s.floor_configuration[floor_number];
            unsigned int unused = 0;

            // printf("%d %d\n", t, fl.type);
            // fflush(stdout);

            if (fl.type == 0) {
                fl.amplitude = in.read_u(fl.amplitude_bits);
                if (fl.amplitude > 0) {
                    int booknumber = in.read_u(ilog(fl.number_of_books));
                    if (booknumber > s.codebook_count)
                        exit(-1);
                    double last = 0;
                    while (1) {
                        vector<double> tmp = s.codebook_configuration[fl.book_list[booknumber]].lookup(in);
                        for (unsigned int i = 0; i < tmp.size(); ++i)
                            tmp[i] += last;
                        last = tmp[tmp.size() - 1];
                        fl.coefficients.insert(fl.coefficients.end(), tmp.begin(), tmp.end());
                        if (fl.coefficients.size() >= (long unsigned int)fl.order)
                            break;
                    }
                    floor0_synthesis(fl, t, n / 2);
                } else
                    unused = 1;
            } else {
                bool nonzero = in.read_u(1);
                if (nonzero == 0) {
                    unused = 1;
                } else {
                    int rr[4] = {256, 128, 86, 64};
                    fl.range = rr[fl.multiplier - 1];
                    fl.Y.clear();
                    fl.Y.resize(fl.values);
                    fl.Y[0] = in.read_u(ilog(fl.range - 1));
                    fl.Y[1] = in.read_u(ilog(fl.range - 1));
                    unsigned int offset = 2;
                    for (int i = 0; i < fl.partitions; ++i) {
                        int class_ = fl.partition_class_list[i];
                        int cdim = fl.class_dimensions[class_];
                        int cbits = fl.class_subclasses[class_];
                        int csub = pow(2, cbits) - 1;
                        int cval = 0;
                        if (cbits > 0)
                            cval = s.codebook_configuration[fl.class_masterbooks[class_]].find(in);
                        for (int j = 0; j < cdim; ++j) {
                            int book = fl.subclass_books[class_][cval & csub];
                            cval >>= cbits;
                            if (book >= 0)
                                fl.Y[j + offset] = s.codebook_configuration[book].find(in);
                            else
                                fl.Y[j + offset] = 0;
                        }
                        offset += cdim;
                    }
                    floor1_synthesis(fl, t, n / 2);
                }
            }

            if (unused == 1) {
                for (int i = 0; i < n / 2; ++i)
                    floor_output[t][i] = 0;
                no_residue[t] = 1;
            }
        }
    }

    void decode_residue(io_buf &in, identification &id, setup &s, modes &mode, mappings &mapping, residues &residue, int ch, bool *do_not_decode_flag) {
        int residue_type = residue.type;
        // decode packet
        int actual_size = n / 2;
        if (residue_type == 2)
            actual_size *= ch;
        int limit_residue_begin = min(residue.begin, actual_size);
        int limit_residue_end = min(residue.end, actual_size);
        int classwords_per_codeword = s.codebook_configuration[residue.classbook].dimensions;
        int n_to_read = limit_residue_end - limit_residue_begin;
        unsigned int partitions_to_read = n_to_read / residue.partition_size + (n_to_read % residue.partition_size != 0);

        r_output.clear();
        r_output_tmp.clear();
        r_output_tmp.resize(ch);
        for (int i = 0; i < ch; ++i)
            r_output_tmp[i].assign(actual_size, 0.0);

        /*
           fprintf(stderr, "%d %d %d %d\n", partitions_to_read, residue.partition_size, actual_size, n);
           fprintf(stderr, "%d %d %d %d\n", limit_residue_begin, limit_residue_end, residue.begin, residue.end);
           fprintf(stderr, "%d %d %d\n", classwords_per_codeword, n_to_read, residue_type);

           for(int i = 0; i < ch; ++i)
           fprintf(stderr, ".%d ", do_not_decode_flag[i]);
           fprintf(stderr, "\n");
         */
        if (in.end_p == 1)
            return;

        if (n_to_read != 0) {

            if (residue_type < 2) {

                int classifications[ch][classwords_per_codeword + partitions_to_read];
                for (int pass = 0; pass < 8; ++pass) {
                    unsigned int partition_count = 0;
                    while (partition_count < partitions_to_read) {

                        if (pass == 0) {
                            for (int j = 0; j < ch; ++j)
                                if (do_not_decode_flag[j] == 0) {
                                    unsigned int temp = s.codebook_configuration[residue.classbook].find(in);
                                    for (int t = classwords_per_codeword - 1; t >= 0; --t) {
                                        classifications[j][t + partition_count] = temp % residue.classification;
                                        temp /= residue.classification;
                                    }
                                }
                        }

                        for (int t = 0; t < classwords_per_codeword && partition_count < partitions_to_read; ++t) {
                            for (int j = 0; j < ch; ++j)
                                if (do_not_decode_flag[j] == 0) {
                                    int vqclass = classifications[j][partition_count];
                                    int vqbook = residue.books[vqclass][pass];
                                    if (vqbook != -1) {
                                        int offset = limit_residue_begin + partition_count * residue.partition_size;
                                        codebooks &code = s.codebook_configuration[vqbook];
                                        if (residue_type == 0) {
                                            int step = residue.partition_size / code.dimensions;
                                            for (int k = 0; k < step; ++k) {
                                                if (in.end_p == 1) {
                                                    partition_count = partitions_to_read, j = ch, pass = 8;
                                                    break;
                                                }
                                                vector<double> entry_temp = code.lookup(in);
                                                for (int l = 0; l < code.dimensions; ++l)
                                                    r_output_tmp[j][offset + k + l * step] += entry_temp[l];
                                            }
                                        } else if (residue_type == 1) {
                                            int k = 0;
                                            do {
                                                if (in.end_p == 1) {
                                                    partition_count = partitions_to_read, j = ch, pass = 8;
                                                    break;
                                                }
                                                vector<double> entry_temp = code.lookup(in);
                                                for (int l = 0; l < code.dimensions; ++l) {
                                                    r_output_tmp[j][offset + k] += entry_temp[l];
                                                    ++k;
                                                }
                                            } while (k < residue.partition_size);
                                        }
                                    }
                                }
                            partition_count++;
                        }
                    }
                }
            } else { // format 2
                vector<double> vv(ch * actual_size, 0);
                bool noskip = 0;
                for (int i = 0; i < ch; ++i)
                    if (do_not_decode_flag[i] == 0)
                        noskip = 1;
                int classifications[classwords_per_codeword + partitions_to_read];
                for (int pass = 0; noskip && pass < 8; ++pass) {
                    unsigned int partition_count = 0;
                    while (partition_count < partitions_to_read) {

                        if (pass == 0) {
                            unsigned int temp = s.codebook_configuration[residue.classbook].find(in);
                            for (int t = classwords_per_codeword - 1; t >= 0; --t) {
                                classifications[t + partition_count] = temp % residue.classification;
                                temp /= residue.classification;
                            }
                        }

                        for (int t = 0; t < classwords_per_codeword && partition_count < partitions_to_read; ++t) {
                            int vqclass = classifications[partition_count];
                            int vqbook = residue.books[vqclass][pass];
                            if (vqbook != -1) {
                                int offset = limit_residue_begin + partition_count * residue.partition_size;
                                codebooks &code = s.codebook_configuration[vqbook];
                                int k = 0;
                                do {
                                    if (in.end_p == 1) {
                                        partition_count = partitions_to_read, pass = 8;
                                        break;
                                    }
                                    vector<double> entry_temp = code.lookup(in);
                                    for (int l = 0; l < code.dimensions; ++l) {
                                        vv[offset + k] += entry_temp[l];
                                        ++k;
                                    }
                                } while (k < residue.partition_size);
                            }
                            partition_count++;
                        }
                    }
                }

                for (int i = 0; i < n_to_read; ++i)
                    for (int j = 0; j < ch; ++j)
                        r_output_tmp[j][i] = vv[i * ch + j];
            }
        }

        /*
           for(int i = 0; i < ch; ++i) {
           r_output_tmp[i].resize(n_to_read);
           }
         */
    }

    void floor0_synthesis(floors &fl, int t, int n) {
        auto foobar = [&](int _i) -> int {
            int ret = (int)(bark(((double)fl.rate) * _i / 2 / n) * fl.bark_map_size / bark(0.5 * fl.rate));
            return ret;
        };
        vector<double> map(n, 0);
        for (int i = 0; i < n - 1; ++i) {
            map[i] = min((int)(fl.bark_map_size - 1), foobar(i));
        }
        map[n - 1] = -1;

        auto getpq = [&](double omega, double &p, double &q) -> void {
            double c = cos(omega);
            p = 1, q = 1;
            if ((fl.order & 1) == 0) {
                for (int j = 0; j <= ((fl.order - 3) / 2); ++j)
                    p *= (fl.coefficients[2 * j + 1] - c) * (fl.coefficients[2 * j + 1 - c] - c);
                for (int j = 0; j <= ((fl.order - 1) / 2); ++j)
                    q *= (fl.coefficients[2 * j] - c) * (fl.coefficients[2 * j] - c);
                p *= (1 - c * c) * 4;
                q /= 4;
            } else {
                for (int j = 0; j <= ((fl.order - 2) / 2); ++j)
                    p *= (fl.coefficients[2 * j + 1] - c) * (fl.coefficients[2 * j + 1 - c] - c);
                for (int j = 0; j <= ((fl.order - 2) / 2); ++j)
                    q *= (fl.coefficients[2 * j] - c) * (fl.coefficients[2 * j] - c);
                p *= (1 - c) * 2;
                q *= (1 + c) * 2;
            }
        };

        auto linear = [&](double p, double q) -> double {
            return exp(.11512925 * (((double)fl.amplitude) * fl.amplitude_offset / (pow(2, fl.amplitude_bits) - 1) / sqrt(p + q) - fl.amplitude_offset));
        };

        for (int i = 0; i < n;) {
            double omega = M_PI * map[i] / fl.bark_map_size;
            double p, q;
            getpq(omega, p, q);
            double linear_floor_value = linear(p, q);
            unsigned int iteration_condition;
            do {
                iteration_condition = map[i];
                floor_output[t][i] = linear_floor_value;
                i++;
            } while (map[i] == iteration_condition);
        }
    }

    void floor1_synthesis(floors &fl, int t, int n) {
        // step 1
        int final_Y[fl.range] = {};
        /*
           printf("%d %d\n", fl.range, fl.values);
           fflush(stdout);
         */
        bool step2_flag[fl.values] = {};
        step2_flag[0] = 1, step2_flag[1] = 1;
        final_Y[0] = fl.Y[0], final_Y[1] = fl.Y[1];
        for (int i = 2; i < fl.values; ++i) {
            int low_neighbor_offset = low_neighbor(fl.X_list, i);
            int high_neighbor_offset = high_neighbor(fl.X_list, i, fl.values);
            int predicted = render_point(
                fl.X_list[low_neighbor_offset], final_Y[low_neighbor_offset],
                fl.X_list[high_neighbor_offset], final_Y[high_neighbor_offset], fl.X_list[i]);
            int val = fl.Y[i];
            int highroom = fl.range - predicted, lowroom = predicted, room;
            if (highroom < lowroom)
                room = highroom * 2;
            else
                room = lowroom * 2;
            if (val != 0) {
                step2_flag[low_neighbor_offset] = 1;
                step2_flag[high_neighbor_offset] = 1;
                step2_flag[i] = 1;
                if (val >= room) {
                    if (highroom > lowroom)
                        final_Y[i] = val - lowroom + predicted;
                    else
                        final_Y[i] = predicted - val + highroom - 1;
                } else {
                    if (val & 1)
                        final_Y[i] = predicted - ((val + 1) / 2);
                    else
                        final_Y[i] = predicted + (val / 2);
                }
            } else {
                step2_flag[i] = 0;
                final_Y[i] = predicted;
            }
        }
        // step 2
        int X_list2[fl.values];
        for (int i = 0; i < fl.values; ++i)
            X_list2[i] = fl.X_list[i];
        floor1_sort(X_list2, final_Y, step2_flag, fl.values);
        /*
           for(int i = 0; i < fl.values; ++i)
           printf("(%d,%d,%d, %d/%d) ", X_list2[i], final_Y[i], step2_flag[i], i, fl.values);
           printf("\n");
           fflush(stdout);
         */
        int hx = 0, lx = 0;
        int ly = final_Y[0] * fl.multiplier, hy = 0;
        int ma_x = 0;
        for (int i = 0; i < fl.values; ++i)
            ma_x = max(ma_x, X_list2[i]);
        vector<int> output(ma_x, 0);
        for (int i = 1; i < fl.values; ++i) {
            /*
               printf("%d/%d\n", i, fl.values);
               fflush(stdout);
             */
            if (step2_flag[i] == 1) {
                hy = final_Y[i] * fl.multiplier;
                hx = X_list2[i];
                // if(lx == hx) exit(-1);
                render_line(lx, ly, hx, hy, output);
                lx = hx, ly = hy;
            }
        }
        if (hx < (int)n)
            render_line(hx, hy, n, hy, output);
        if (hx > (int)n)
            output.resize(n);
        for (int i = 0; i < n; ++i) {
            floor_output[t][i] = floor1_inverse_dB_table[output[i]];
        }
        output.clear();
    }

    void decode(io_buf &in, identification &id, setup &s) {
        // packet type, mode and window decode
        /*
           printf("window\n");
           fflush(stdout);
         */
        mode_number = in.read_u(ilog(s.mode_count - 1));
        auto &mode = s.mode_configuration[mode_number];
        mapping_number = mode.mapping;
        auto &mapping = s.mapping_configuration[mapping_number];

        decode_window(in, id, s, mode, mapping);

        /*
           printf("floor_curve\n");
           fflush(stdout);
         */
        decode_floor_curve(in, id, s, mode, mapping);

        // nonzero vector propagate
        /*
           printf("propagate\n");
           fflush(stdout);
         */
        vector<vector<double>> residue_output(id.audio_channels);
        for (int i = 0; i < mapping.submap; ++i) {
            int ch = 0;
            bool do_not_decode_flag[id.audio_channels] = {};
            for (unsigned int j = 0; j < id.audio_channels; ++j) {
                if (mapping.mux[j] == i) {
                    if (no_residue[j] == 1)
                        do_not_decode_flag[ch] = 1;
                    else
                        do_not_decode_flag[ch] = 0;
                    ch++;
                }
            }
            int residue_number = mapping.submap_residue[i];
            auto &residue = s.residue_configuration[residue_number];
            decode_residue(in, id, s, mode, mapping, residue, ch, do_not_decode_flag);
            ch = 0;
            for (unsigned int j = 0; j < id.audio_channels; ++j) {
                if (mapping.mux[j] == i) {
                    residue_output[j].assign(r_output_tmp[ch].begin(), r_output_tmp[ch].end());
                    /*
                       for(int k = 0; k < residue_output[j].size(); ++k)
                       residue_output[j][k] = max(min(residue_output[j][k], 1.0), -1.0);
                     */
                    ch++;
                }
            }
        }

        // inverse coupling
        for (int i = mapping.coupling_step - 1; i >= 0; --i) {
            vector<double> &magnitude = residue_output[mapping.magnitude[i]];
            vector<double> &angle = residue_output[mapping.angle[i]];
            for (int j = 0; j < n / 2; ++j) {
                double M = magnitude[j], A = angle[j], new_M, new_A;
                if (M > 0 && A > 0)
                    new_M = M, new_A = M - A;
                else if (M > 0 && A <= 0)
                    new_A = M, new_M = M + A;
                else if (M <= 0 && A > 0)
                    new_M = M, new_A = M + A;
                else
                    new_A = M, new_M = M - A;
                magnitude[j] = new_M, angle[j] = new_A;
            }
        }

        // dot product
        /*
           printf("dot product\n");
           fflush(stdout);
         */
        for (unsigned int i = 0; i < id.audio_channels; ++i) {
            for (int j = 0; j < n / 2; ++j) {
                // printf("%.2lf,", floor_output[i][j]);
                floor_output[i][j] *= round(residue_output[i][j]);
                // fprintf(stderr, "%.5f ", round(residue_output[i][j]));
            }
            // fprintf(stderr, "\n");
        }

        // inverse MDCT
        /*
           printf("inverse MDCT\n");
           fflush(stdout);
         */
        for (int i = 0; i < id.audio_channels; ++i) {
            Y[i].clear();
            Y[i].assign(n, 0);
            easy_IMDCT(floor_output[i], Y[i], window, n);
        }
        /*
           for(int i = 0; i < id.audio_channels; ++i) {
           for(int j =0; j < n; ++j)
           fprintf(stderr, "%.5lf ", Y[i][j]);
           fprintf(stderr, "\n");
           }
         */

        // overlap_add
        /*
           printf("overlap_add\n");
           fflush(stdout);
         */
        int m, m2;
        m = (previous_window_flag == 0) ? id.blocksize_0 : id.blocksize_1;
        m2 = m * 3 / 4;
        int l2 = left_n / 2;
        int left_l = m2 - l2;
        /*
           printf("%d %d\n", previous_window_flag, next_window_flag);
           printf("%d %d %d %d\n", n, m, m2, l2);
           printf("%d %d %d %d\n", left_l, left_window_start, left_window_end, left_n);
           printf("%d\n", m / 4 + n / 4);
         */
        for (int i = 0; i < id.audio_channels; ++i) {
            Y_ret[i].clear();
            for (int j = m / 2; j < left_l; ++j)
                Y_ret[i].push_back(previous_y[i][j]);
            for (int j = 0; j < left_n; ++j)
                Y_ret[i].push_back(previous_y[i][j + left_l] + Y[i][j + left_window_start]);
            for (int j = left_window_end; j < n / 2; ++j)
                Y_ret[i].push_back(Y[i][j]);
        }

        /*
           printf("assign result\n");
           fflush(stdout);
         */
        for (int i = 0; i < id.audio_channels; ++i) {
            previous_y[i].clear();
            previous_y[i].assign(Y[i].begin(), Y[i].end());
        }
    }
};

#endif
