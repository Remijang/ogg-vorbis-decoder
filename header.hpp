#ifndef HEADER_HPP
#define HEADER_HPP

#include <vector>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <set>
#include "io.hpp"
#include "util.hpp"
#include "huffman.hpp"
#include "mdct.hpp"

using namespace std;

struct identification {
	unsigned int  verbis_version;
	unsigned char audio_channels;
	unsigned int  audio_sample_rate;
	int			  bitrate_maximum;
	int			  bitrate_nominal;
	int			  bitrate_minimum;
	unsigned int  blocksize_0;
	unsigned int  blocksize_1;
	bool		  framing_flag;

	void init (io_buf &in) {
		verbis_version    = in.read_u(32);
		audio_channels    = in.read_u(8);
		audio_sample_rate = in.read_u(32); 
		bitrate_maximum   = in.read_s(32);
		bitrate_nominal   = in.read_s(32);
		bitrate_minimum   = in.read_s(32);
		blocksize_0       = pow(2, (int)in.read_u(4));
		blocksize_1       = pow(2, (int)in.read_u(4));
		framing_flag      = in.read_u(1);
	}
};

struct comment {
	unsigned int vender_length;
	vector<unsigned char> vender_string;
	unsigned int user_comment_list_length;
	vector<vector<unsigned char>> user_comment_list;
	bool framing_bit;

	vector<unsigned char> get_utf_string(io_buf &in, unsigned int n) {
		vector<unsigned char> ret;
		unsigned int i = 0;
		for(i = 0; i < n; ) {
			vector<unsigned char> utf = in.read_utf();
			ret.insert(ret.end(), utf.begin(), utf.end());
			i += utf.size();
		}
		if(i != n) exit(-1);
		return ret;
	}

	void init (io_buf &in) {
		vender_length = in.read_u(32);
		vender_string = get_utf_string(in, vender_length);
		user_comment_list_length = in.read_u(32);
		for(unsigned int i = 0; i < user_comment_list_length; ++i) {
			unsigned int length = in.read_u(32);
			user_comment_list.push_back(get_utf_string(in, length));
		}
		framing_bit = in.read_u(1);
	}

};

struct codebooks {
	int dimensions;
	int entries;
	int* length;
	long long* codeword;
	bool ordered;
	int lookup_type;
	int* multiplicands;
	double minimum_value;
	double delta_value;
	int value_bits;
	bool sequence_p;
	int lookup_values;

	~codebooks() {
		delete[] length;
		delete[] codeword;
		delete[] multiplicands;
	}

	unsigned int find(io_buf &in) {
		int ok = 0;
		long long s = 0;
		int current_length = 0;
		while(!ok) {
			s <<= 1, s += in.read_u(1), current_length++;
			for(int i = 0; i < entries; ++i) { 
				if(current_length == length[i] && s == codeword[i]) return i;
			}
			if(current_length > 64){
				printf("not found\n");
				exit(-1);
			}
		}
		return 0;
	}

	vector<double> lookup(io_buf &in) {
		unsigned int lookup_offset = find(in);
		vector<double> value_vector(dimensions, 0);
		if(lookup_type == 1) {
			double last = 0;
			int index_divisor = 1;
			for(int i = 0; i < dimensions; ++i) {
				unsigned int multiplicand_offset = (lookup_offset / index_divisor) % lookup_values;
				value_vector[i] = multiplicands[multiplicand_offset] * delta_value + minimum_value + last;
				if(sequence_p == 1) last = value_vector[i];
				index_divisor *= lookup_values;
			}
		}
		if(lookup_type == 2) {
			double last = 0;
			int multiplicand_offset = lookup_offset * dimensions;
			for(int i = 0; i < dimensions; ++i) {
				value_vector[i] = multiplicands[multiplicand_offset] * delta_value + minimum_value + last;
				if(sequence_p == 1) last = value_vector[i];
				multiplicand_offset++;
			}
		}
		return value_vector;
	}

	void decode(io_buf &in) {
		unsigned int tmp = in.read_u(24);
		if(tmp != 0x564342) exit(-1);
		dimensions = in.read_u(16);
		entries    = in.read_u(24);
		ordered    = in.read_u(1);
		length = new int[entries];
		codeword = new long long[entries];
		if(ordered == 0) {
			bool sparse = in.read_u(1);
			for(int i = 0; i < entries; ++i) {
				if(sparse == 1) {
					bool flag = in.read_u(1);
					if(flag == 1) length[i] = in.read_u(5) + 1;
					else length[i] = 0;
				}
				else {
					length[i] = in.read_u(5) + 1;
				}
			}
		}
		else {
			int current_entry = 0;
			int current_length = in.read_u(5) + 1;
			while(1) {
				unsigned int number = in.read_u(ilog(entries - current_entry));
				for(unsigned int i = 0; i < number; ++i)
					length[i + current_entry] = current_length;
				current_entry += number;
				current_length++;
				if(current_entry > entries) exit(-1);
				else if(current_entry == entries) break;
			}
		}
		huffman(entries, length, codeword);
		/*
		for(int i = 0; i < entries; ++i) {
			long long u = codeword[i];
			for(int t = length[i] - 1; t >= 0; --t) {
				printf("%d", (u >> t) & 1);
			}
			printf(" %d\n", length[i]);
		}*/

		lookup_type = in.read_u(4);
		if(lookup_type == 0) {
			return;
		}
		minimum_value = float32_unpack(in.read_u(32));
		delta_value   = float32_unpack(in.read_u(32));
		value_bits = in.read_u(4) + 1;
		sequence_p = in.read_u(1);
		if(lookup_type == 1) lookup_values = lookup1_values(entries, dimensions);
		else if(lookup_type == 2) lookup_values = entries * dimensions;
		else exit(-1);

		multiplicands = new int[lookup_values];
		for(int i = 0; i < lookup_values; ++i)
			multiplicands[i] = in.read_u(value_bits);

		return;
	}
};

struct floors {
	int type;
	// type 0 //
	int order;
	int rate;
	int bark_map_size;
	int amplitude_bits;
	int amplitude_offset;
	int number_of_books;
	int* book_list;
	vector<double> coefficients;
	unsigned long long amplitude;

	// type 1 //
	int partitions;
	int maximum_class;
	int* partition_class_list;
	int* class_dimensions;
	int* class_subclasses;
	int* class_masterbooks;
	int** subclass_books;
	int multiplier;
	int rangebits;
	int* X_list;
	int values;
	int range;
	vector<int> Y;

	~floors() {
		if(type== 0) delete[] book_list;
		else {
			delete[] partition_class_list;
			delete[] class_dimensions;
			delete[] class_subclasses;
			delete[] class_masterbooks;
			for(int i = 0; i < maximum_class; ++i)
				delete[] subclass_books[i];
			delete[] subclass_books;
			delete[] X_list;
		}
	}

	void header_decode (io_buf &in, unsigned int _type, int co) {
		type = _type;
		if(type == 0) {
			order = in.read_u(8);
			rate = in.read_u(16);
			bark_map_size = in.read_u(16);
			amplitude_bits = in.read_u(6);
			amplitude_offset = in.read_u(8);
			number_of_books = in.read_u(4) + 1;
			book_list = new int[number_of_books];
			for(int i = 0; i < number_of_books; ++i)
				book_list[i] = in.read_u(8);
		}
		else if(type == 1) {
			partitions = in.read_u(5);
			maximum_class = -1;
			partition_class_list = new int[partitions];
			for(int i = 0; i < partitions; ++i) {
				int tmp = in.read_u(4);
				if(maximum_class < (int) tmp) maximum_class = tmp;
				partition_class_list[i] = tmp;
			}
			maximum_class++;
			class_dimensions  = new int[maximum_class];
			class_subclasses  = new int[maximum_class];
			class_masterbooks = new int[maximum_class];
			subclass_books = new int*[maximum_class];
			for(int i = 0; i < maximum_class; ++i) {
				class_dimensions[i] = in.read_u(3) + 1;
				class_subclasses[i] = in.read_u(2);
				if(class_subclasses[i] != 0) class_masterbooks[i] = in.read_u(8);
				int p2 = pow(2, class_subclasses[i]);
				subclass_books[i] = new int[p2];
				for(int j = 0; j < p2; ++j) {
					subclass_books[i][j] = ((int) in.read_u(8)) - 1;
					if(subclass_books[i][j] >= co) exit(-1);
				}
			}
			multiplier = in.read_u(2) + 1;
			rangebits = in.read_u(4);
			int total = 0;
			for(int i = 0; i < partitions; ++i)
				total += class_dimensions[partition_class_list[i]];
			X_list = new int[total + 2];
			X_list[0] = 0;
			X_list[1] = pow(2, rangebits);
			values = 2;
			for(int i = 0; i < partitions; ++i) {
				int current_class_number = partition_class_list[i];
				for(int j = 0; j < class_dimensions[current_class_number]; ++j) {
					X_list[values] = in.read_u(rangebits);
					values++;
				}
			}
		}
		return;
	}
};

struct residues {
	int type;
	int begin;
	int end;
	int partition_size;
	int classification;
	int classbook;
	int* cascade;
	int** books;

	~residues() {
		delete[] cascade;
		for(int i = 0; i < classification; ++i)
			delete[] books[i];
		delete[] books;
	}

	void header_decode(io_buf &in, int _type) {
		type = _type;
		begin = in.read_u(24);
		end = in.read_u(24);
		partition_size = in.read_u(24) + 1;
		classification = in.read_u(6) + 1;
		classbook = in.read_u(8);
		cascade = new int[classification];
		for(int i = 0; i < classification; ++i) {
			int high_bits = 0;
			int low_bits = in.read_u(3);
			bool bitflag = in.read_u(1);
			if(bitflag == 1) high_bits = in.read_u(5);
			cascade[i] = high_bits * 8 + low_bits;
		}
		books = new int*[classification];
		for(int i = 0; i < classification; ++i) {
			books[i] = new int[8];
			for(int j = 0; j < 8; ++j) {
				if((cascade[i] >> j) & 1) books[i][j] = in.read_u(8);
				else books[i][j] = -1; // unused
			}
		}
		return;
	}
};

struct mappings {
	bool flag;
	int coupling_step;
	int submap;
	int* magnitude;
	int* angle;
	int* mux;
	int* submap_floor;
	int* submap_residue;

	~mappings() {
		if(flag == 1) {
			delete[] magnitude;
			delete[] angle;
		}
		delete[] mux;
		delete[] submap_floor;
		delete[] submap_residue;
	}

	void header_decode (io_buf &in, identification &id, int fl, int re) {
		flag = in.read_u(1);
		if(flag == 1) submap = in.read_u(4) + 1;
		else submap = 1;
		flag = in.read_u(1);
		if(flag == 1) {
			coupling_step = in.read_u(8) + 1;
			magnitude = new int[coupling_step];
			angle = new int[coupling_step];
			for(int i = 0; i < coupling_step; ++i) {
				magnitude[i] = in.read_u(ilog(id.audio_channels - 1));
				angle[i] = in.read_u(ilog(id.audio_channels - 1));
				if(angle[i] == magnitude[i] 
						|| magnitude[i] > id.audio_channels - 1 
						|| angle[i] > id.audio_channels - 1) exit(-1);
			}
		}
		else coupling_step = 0;
		int rev = in.read_u(2);
		if(rev != 0) exit(-1);
		if(submap > 1) {
			mux = new int[id.audio_channels];
			for(int i = 0; i < id.audio_channels; ++i) {
				mux[i] = in.read_u(4);
				if(mux[i] > submap - 1) exit(-1);
			}
		}
		else {
			mux = new int[id.audio_channels];
			for(int i = 0; i < id.audio_channels; ++i) {
				mux[i] = 0;
			}
		}
		submap_floor = new int[submap];
		submap_residue = new int [submap];
		for(int i = 0; i < submap; ++i) {
			in.read_u(8); // unused
			submap_floor[i] = in.read_u(8);
			if(submap_floor[i] >= fl) exit(-1);
			submap_residue[i] = in.read_u(8);
			if(submap_residue[i] >= re) exit(-1);
		}
		return;
	}
};

struct modes {
	bool blockflag;
	int windowtype;
	int transformtype;
	int mapping;

	void header_decode (io_buf &in, int ma) {
		blockflag = in.read_u(1);
		windowtype = in.read_u(16);
		transformtype = in.read_u(16);
		mapping = in.read_u(8);
		if(windowtype != 0
				|| transformtype != 0 
				|| mapping > ma) exit(-1);
		return;
	}
};

struct setup {
	int codebook_count;
	codebooks* codebook_configuration;
	int time_count;
	int* time;
	int floor_count;
	int* floor_type;
	floors* floor_configuration;
	int residue_count;
	int* residue_type;
	residues* residue_configuration;
	int mapping_count;
	mappings* mapping_configuration;
	int mode_count;
	modes* mode_configuration;

	~setup() {
		delete[] codebook_configuration;
		delete[] time;
		delete[] floor_type;
		delete[] floor_configuration;
		delete[] residue_type;
		delete[] residue_configuration;
		delete[] mapping_configuration;
		delete[] mode_configuration;
	}

	void init(io_buf &in, identification &id) {
		codebook_count = in.read_u(8) + 1;
		codebook_configuration = new codebooks[codebook_count];
		for(int i = 0; i < codebook_count; ++i) {
			codebook_configuration[i].decode(in);
		}

		time_count = in.read_u(6) + 1;
		time = new int[time_count];
		for(int i = 0; i < time_count; ++i) {
			time[i] = in.read_u(16);
			if(time[i] != 0) exit(-1);
		}

		floor_count = in.read_u(6) + 1;
		floor_type = new int[floor_count];
		floor_configuration = new floors[floor_count];
		for(int i = 0; i < floor_count; ++i) {
			floor_type[i] = in.read_u(16);
			if(floor_type[i] > 1) exit(-1);
			floor_configuration[i].header_decode(in, floor_type[i], codebook_count);
		}
		
		residue_count = in.read_u(6) + 1;
		residue_type = new int[residue_count];
		residue_configuration = new residues[residue_count];
		for(int i = 0; i < residue_count; ++i) {
			residue_type[i] = in.read_u(16);
			if(residue_type[i] > 2) exit(-1);
			residue_configuration[i].header_decode(in, residue_type[i]);
		}

		mapping_count = in.read_u(6) + 1;
		mapping_configuration = new mappings[mapping_count];
		for(int i = 0; i < mapping_count; ++i) {
			int mapping_type = in.read_u(16);
			if(mapping_type != 0) exit(-1);
			mapping_configuration[i].header_decode(in, id, floor_count, residue_count);
		}

		mode_count = in.read_u(6) + 1;
		mode_configuration = new modes[mode_count];
		for(int i = 0 ; i < mode_count; ++i) {
			mode_configuration[i].header_decode(in, mapping_count);
		}
		bool framing = in.read_u(1);
		if(framing == 0) exit(-1);
	}

};

struct packet {
	//unsigned int type;
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
		for(int i = 0; i < id.audio_channels; ++i)
			previous_y[i].assign(n, 0);
	}

	void decode_window(io_buf &in, identification &id, setup &s, modes &mode, mappings& mapping) {
		unsigned int blockflag = mode.blockflag;
		previous_window_flag = next_window_flag = 0;
		n = (blockflag == 0) ? (int) id.blocksize_0 : (int) id.blocksize_1;
		if(blockflag == 1) {
			previous_window_flag = in.read_u(1);
			next_window_flag = in.read_u(1);
		}

		window.clear();
		window.resize(n);
		window_center = n / 2;
		if(blockflag == 1 && previous_window_flag == 0) {
			left_window_start = n / 4 - id.blocksize_0 / 4;
			left_window_end = n / 4 + id.blocksize_0 / 4;
			left_n = id.blocksize_0 / 2;
		}
		else {
			left_window_start = 0;
			left_window_end = window_center;
			left_n = n / 2;
		}
		if(blockflag == 1 && next_window_flag == 0) {
			right_window_start = n * 3 / 4 - id.blocksize_0 / 4;
			right_window_end = n * 3 / 4 + id.blocksize_0 / 4;
			right_n = id.blocksize_0 / 2;
		}
		else {
			right_window_start = window_center;
			right_window_end = n;
			right_n = n / 2;
		}
		for(int i = 0; i < left_window_start; ++i) window[i] = 0;
		for(int i = left_window_start; i < left_window_end; ++i)
			window[i] = sin(M_PI / 2 * pow(sin((((double)i) - ((double)left_window_start) + 0.5) / left_n * ( M_PI / 2)), 2));
		for(int i = right_window_start; i < right_window_end; ++i)
			window[i] = sin(M_PI / 2 * pow(sin((((double)i) - ((double)right_window_start) + 0.5) / right_n * ( M_PI / 2) + M_PI / 2), 2));
		for(int i = right_window_end; i < n; ++i) window[i] = 0;
	}

	void decode_floor_curve(io_buf &in, identification &id, setup &s, modes &mode, mappings &mapping) {
		// floor curve decode
		// vector<vector<double>> floor_output;
		floor_output.clear();
		floor_output.resize(id.audio_channels);
		for(unsigned int t = 0; t < id.audio_channels; ++t)
			floor_output[t].assign(n / 2, 0);

		no_residue.clear();
		no_residue.resize(id.audio_channels);
		for(int i = 0; i < id.audio_channels; ++i) no_residue[i] = 0;

		for(unsigned int t = 0; t < id.audio_channels; ++t) {
			/*
			unsigned int submap_number;
			if(mapping.submap > 1) submap_number = mapping.mux[t];
			else submap_number = 0;
			*/
			unsigned int submap_number = mapping.mux[t];
			unsigned int floor_number = mapping.submap_floor[submap_number];
			auto &fl = s.floor_configuration[floor_number];
			unsigned int unused = 0;

			//printf("%d %d\n", t, fl.type);
			//fflush(stdout);

			if(fl.type == 0) {
				fl.amplitude = in.read_u(fl.amplitude_bits);
				if(fl.amplitude > 0) {
					int booknumber = in.read_u(ilog(fl.number_of_books));
					if(booknumber > s.codebook_count) exit(-1);
					double last = 0;
					while(1) {
						vector<double> tmp = s.codebook_configuration[fl.book_list[booknumber]].lookup(in);
						for(unsigned int i = 0; i < tmp.size(); ++i) tmp[i] += last;
						last = tmp[tmp.size() - 1];
						fl.coefficients.insert(fl.coefficients.end(), tmp.begin(), tmp.end());
						if(fl.coefficients.size() >= (long unsigned int) fl.order) break;
					}
					floor0_synthesize(fl, t, n / 2);
				}
				else unused = 1;
			}
			else {
				bool nonzero = in.read_u(1);
				if(nonzero == 0) {
					unused = 1;
				}
				else {
					int rr[4] = {256, 128, 86, 64};
					fl.range = rr[fl.multiplier - 1];
					fl.Y.clear();
					fl.Y.resize(fl.values);
					fl.Y[0] = in.read_u(ilog(fl.range - 1));
					fl.Y[1] = in.read_u(ilog(fl.range - 1));
					unsigned int offset = 2;
					for(int i = 0; i < fl.partitions; ++i) {
						int class_ = fl.partition_class_list[i];
						int cdim   = fl.class_dimensions[class_];
						int cbits  = fl.class_subclasses[class_];
						int csub   = pow(2, cbits) - 1;
						int cval   = 0;
						if(cbits > 0) cval = s.codebook_configuration[fl.class_masterbooks[class_]].find(in);
						for(int j = 0; j < cdim; ++j) {
							int book = fl.subclass_books[class_][cval & csub];
							cval >>= cbits;
							if(book >= 0) fl.Y[j + offset] = s.codebook_configuration[book].find(in);
							else fl.Y[j + offset] = 0;
						}
						offset += cdim;
					}
					floor1_synthesize(fl, t, n / 2);
				}
			}

			if(unused == 1) {
				for(int i = 0; i < n / 2; ++i) floor_output[t][i] = 0;
				no_residue[t] = 1;
			}
		}
	}

	void decode_residue(io_buf &in, identification &id, setup &s, modes &mode, mappings &mapping, residues &residue, int ch, bool *do_not_decode_flag) {
		int residue_type = residue.type;
		// decode packet
		int actual_size = n / 2;
		if(residue_type == 2) actual_size *= ch;
		int limit_residue_begin = min(residue.begin, actual_size);
		int limit_residue_end = min(residue.end, actual_size);
		int classwords_per_codeword = s.codebook_configuration[residue.classbook].dimensions;
		int n_to_read = limit_residue_end - limit_residue_begin;
		unsigned int partitions_to_read = n_to_read / residue.partition_size + (n_to_read % residue.partition_size != 0);

		r_output.clear();
		r_output_tmp.resize(ch);
		for(int i = 0; i < ch; ++i)
			r_output_tmp[i].resize(actual_size, 0);

		/*
		printf("%d %d %d %d\n", partitions_to_read, residue.partition_size, actual_size, n);
		printf("%d %d %d %d\n", limit_residue_begin, limit_residue_end, residue.begin, residue.end);
		printf("%d %d %d\n", classwords_per_codeword, n_to_read, residue_type);

		for(int i = 0; i < ch; ++i)
			printf(".%d ", do_not_decode_flag[i]);
		printf("\n");
		*/
		if(in.end_p == 1) return;

		if(n_to_read != 0) {

			if(residue_type < 2) {

				int classifications[ch][classwords_per_codeword + partitions_to_read];
				for(int pass = 0; pass < 8; ++pass) {
					unsigned int partition_count = 0;
					while(partition_count < partitions_to_read) {

						if(pass == 0) {
							for(int j = 0; j < ch; ++j) if(do_not_decode_flag[j] == 0) {
								unsigned int temp = s.codebook_configuration[residue.classbook].find(in);
								for(int t = classwords_per_codeword - 1; t >= 0; --t) {
									classifications[j][t + partition_count] = temp % residue.classification;
									temp /= residue.classification;
								}
							}
						}

						for(int t = 0; t < classwords_per_codeword && partition_count < partitions_to_read; ++t) {
							for(int j = 0; j < ch; ++j) if(do_not_decode_flag[j] == 0) {
								int vqclass = classifications[j][partition_count];
								int vqbook = residue.books[vqclass][pass];
								if(vqbook != -1) {
									int offset = limit_residue_begin + partition_count * residue.partition_size;
									codebooks &code = s.codebook_configuration[vqbook];
									if(residue_type == 0) {
										int step = residue.partition_size / code.dimensions;
										for(int k = 0; k < step; ++k) {
											if(in.end_p == 1) {
												partition_count = partitions_to_read, j = ch, pass = 8;
												break;
											}
											vector<double> entry_temp = code.lookup(in);
											for(int l = 0; l < code.dimensions; ++l) r_output_tmp[j][offset + k + l * step] += entry_temp[l];
										}
									}
									else if(residue_type == 1) {
										int k = 0;
										do {
											if(in.end_p == 1) {
												partition_count = partitions_to_read, j = ch, pass = 8;
												break;
											}
											vector<double> entry_temp = code.lookup(in);
											for(int l = 0; l < code.dimensions; ++l) {
												r_output_tmp[j][offset + k] += entry_temp[l];
												++k;
											}
										} while(k < residue.partition_size);
									}
								}
							}
							partition_count++;
						}

					}
				}
			}
			else { // format 2
				vector<double> vv(ch * actual_size, 0);
				bool noskip = 0;
				for(int i = 0; i < ch; ++i) if(do_not_decode_flag[i] == 0) noskip = 1;
				int classifications[ch][classwords_per_codeword + partitions_to_read];
				for(int pass = 0; noskip && pass < 8; ++pass) {
					unsigned int partition_count = 0;
					while(partition_count < partitions_to_read) {
						
						if(pass == 0) {
							for(int j = 0; j < ch; ++j) {
								unsigned int temp = s.codebook_configuration[residue.classbook].find(in);
								for(int t = classwords_per_codeword - 1; t >= 0; --t) {
									classifications[j][t + partition_count] = temp % residue.classification;
									temp /= residue.classification;
								}
							}
						}

						for(int t = 0; t < classwords_per_codeword && partition_count < partitions_to_read; ++t) {
							for(int j = 0; j < ch; ++j) {
								int vqclass = classifications[j][partition_count];
								int vqbook = residue.books[vqclass][pass];
								if(vqbook != -1) {
									int offset = limit_residue_begin + partition_count * residue.partition_size;
									codebooks &code = s.codebook_configuration[vqbook];
									int k = 0;
									do {
										if(in.end_p == 1) {
											partition_count = partitions_to_read, j = ch, pass = 8;
											break;
										}
										vector<double> entry_temp = code.lookup(in);
										for(int l = 0; l < code.dimensions; ++l) {
											vv[j * residue.partition_size + offset + k] += entry_temp[l];
											++k;
										}
									} while(k < residue.partition_size);
								}
							}
							partition_count++;
						}

					}
				}

				for(int i = 0; i < n_to_read; ++i)
					for(int j = 0; j < ch; ++j)
						r_output_tmp[j][i] = vv[i * ch + j];

			}

		}

		/*
		for(int i = 0; i < ch; ++i) {
			r_output_tmp[i].resize(n_to_read);
		}
		*/
	}

	void floor0_synthesize(floors &fl, int t, int n) {
		auto foobar = [&] (int _i) -> int {
			int ret = (int) (bark(((double) fl.rate) * _i / 2 / n) * fl.bark_map_size / bark(0.5 * fl.rate));
			return ret;
		};
		vector<double> map(n, 0);
		for(int i = 0; i < n - 1; ++i) {
			map[i] = min((int) (fl.bark_map_size - 1), foobar(i));
		}
		map[n - 1] = -1;

		auto getpq = [&] (double omega, double &p, double &q) -> void {
			double c = cos(omega);
			p = 1, q = 1;
			if((fl.order & 1) == 0) {
				for(int j = 0; j <= ((fl.order - 3) / 2); ++j)
					p *= (fl.coefficients[2 * j + 1] - c) * (fl.coefficients[2 * j + 1 - c] - c);
				for(int j = 0; j <= ((fl.order - 1) / 2); ++j)
					q *= (fl.coefficients[2 * j] - c) * (fl.coefficients[2 * j] - c);
				p *= (1 - c * c) * 4;
				q /= 4;
			}
			else {
				for(int j = 0; j <= ((fl.order - 2) / 2); ++j)
					p *= (fl.coefficients[2 * j + 1] - c) * (fl.coefficients[2 * j + 1 - c] - c);
				for(int j = 0; j <= ((fl.order - 2) / 2); ++j)
					q *= (fl.coefficients[2 * j] - c) * (fl.coefficients[2 * j] - c);
				p *= (1 - c) * 2;
				q *= (1 + c) * 2;
			}
		};

		auto linear = [&] (double p, double q) -> double {
			return exp(.11512925 * (((double) fl.amplitude) * fl.amplitude_offset / (pow(2, fl.amplitude_bits) - 1) / sqrt(p + q) - fl.amplitude_offset));
		};

		for(int i = 0; i < n; ) {
			double omega = M_PI * map[i] / fl.bark_map_size;
			double p, q;
			getpq(omega, p, q);
			double linear_floor_value = linear(p, q);
			unsigned int iteration_condition;
			do {
				iteration_condition = map[i];
				floor_output[t][i] = linear_floor_value;
				i++;
			} while(map[i] == iteration_condition);
		}
	}

	void floor1_synthesize(floors &fl, int t, int n) {
		// step 1
		int final_Y[fl.range] = {};
		/*
		printf("%d %d\n", fl.range, fl.values);
		fflush(stdout);
		*/
		bool step2_flag[fl.values] = {};
		step2_flag[0] = 1, step2_flag[1] = 1;
		final_Y[0] = fl.Y[0], final_Y[1] = fl.Y[1];
		for(int i = 2; i < fl.values; ++i) {
			int low_neighbor_offset = low_neighbor(fl.X_list, i);
			int high_neighbor_offset = high_neighbor(fl.X_list, i, fl.values);
			int predicted = render_point(
					fl.X_list[low_neighbor_offset], final_Y[low_neighbor_offset],
					fl.X_list[high_neighbor_offset], final_Y[high_neighbor_offset], fl.X_list[i]);
			int val = fl.Y[i];
			int highroom = fl.range - predicted, lowroom = predicted, room;
			if(highroom < lowroom) room = highroom * 2;
			else room = lowroom * 2;
			if(val != 0) {
				step2_flag[low_neighbor_offset] = 1;
				step2_flag[high_neighbor_offset] = 1;
				step2_flag[i] = 1;
				if(val >= room){
					if(highroom > lowroom) final_Y[i] = val - lowroom + predicted;
					else final_Y[i] = predicted - val + highroom - 1;
				}
				else {
					if(val & 1) final_Y[i] = predicted - ((val + 1) / 2);
					else final_Y[i] = predicted + (val / 2);
				}
			}
			else {
				step2_flag[i] = 0;
				final_Y[i] = predicted;
			}
		}
		// step 2
		int X_list2[fl.values];
		for(int i = 0; i < fl.values; ++i) X_list2[i] = fl.X_list[i];
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
		for(int i = 0; i < fl.values; ++i) ma_x = max(ma_x, X_list2[i]);
		vector<int> output(ma_x, 0);
		for(int i = 1; i < fl.values; ++i) {
			/*
			printf("%d/%d\n", i, fl.values);
			fflush(stdout);
			*/
			if(step2_flag[i] == 1) {
				hy = final_Y[i] * fl.multiplier;
				hx = X_list2[i];
				//if(lx == hx) exit(-1);
				render_line(lx, ly, hx, hy, output);
				lx = hx, ly = hy;
			}
		}
		if(hx < (int) n) render_line(hx, hy, n, hy, output);
		if(hx > (int) n) output.resize(n);
		for(int i = 0; i < n; ++i) {
			floor_output[t][i] = floor1_inverse_dB_table[output[i]];
		}
		output.clear();

	}

	void decode(io_buf &in, identification &id, setup &s) {
		//packet type, mode and window decode
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
		for(int i = 0; i < mapping.submap; ++i) {
			int ch = 0;
			bool do_not_decode_flag[id.audio_channels] = {};
			for(unsigned int j = 0; j < id.audio_channels; ++j) {
				if(mapping.mux[j] == i) {
					if(no_residue[j] == 1) do_not_decode_flag[ch] = 1;
					else do_not_decode_flag[ch] = 0;
					ch++;
				}
			}
			int residue_number = mapping.submap_residue[i];
			auto &residue = s.residue_configuration[residue_number];
			decode_residue(in, id, s, mode, mapping, residue, ch, do_not_decode_flag);
			ch = 0;
			for(unsigned int j = 0; j < id.audio_channels; ++j) {
				if(mapping.mux[j] == i) {
					residue_output[j].assign(r_output_tmp[ch].begin(), r_output_tmp[ch].end());
					for(int k = 0; k < residue_output[j].size(); ++k)
						residue_output[j][k] = max(min(residue_output[j][k], 1.0), -1.0);
					ch++;
				}
			}
		}

		// inverse coupling
		vector<int> magnitude;
		vector<int> angle;
		if(mapping.coupling_step) {
			for(int i = 0; i < mapping.coupling_step; ++i){
				magnitude.push_back(mapping.magnitude[i]);
				angle.push_back(mapping.angle[i]);
			}
		}
		for(int i = mapping.coupling_step - 1; i >= 0; --i) {
			int M = magnitude[i], A = angle[i], new_M, new_A;
			if     (M > 0 && A > 0)  new_M = M, new_A = M - A;
			else if(M > 0 && A <= 0) new_A = M, new_M = M + A;
			else if(M <= 0 && A > 0) new_M = M, new_A = M + A;
			else					 new_A = M, new_M = M - A;
			magnitude[i] = new_M, angle[i] = new_A;
		}

		// dot product
		/*
		printf("dot product\n");
		fflush(stdout);
		*/
		for(unsigned int i = 0; i < id.audio_channels; ++i) {
			for(int j = 0; j < n / 2; ++j) {
				//printf("(%.2lf,", floor_output[i][j]);
				floor_output[i][j] *= residue_output[i][j];
				//printf("%lf) ", residue_output[i][j]);
			}
			//printf("\n");
		}
		//printf("\n");

		// inverse MDCT
		/*
		printf("inverse MDCT\n");
		fflush(stdout);
		*/
		for(int i = 0; i < id.audio_channels; ++i) {
			Y[i].clear();
			Y[i].assign(n, 0);
			easy_IMDCT(floor_output[i], Y[i], window, n);
		}
		/*
		for(int i  =0; i < n; ++i)
			printf("%lf ", Y[0][i]);
		printf("\n");
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
		for(int i = 0; i < id.audio_channels; ++i) {
			Y_ret[i].clear();
			for(int j = m / 2; j < left_l; ++j)
				Y_ret[i].push_back(previous_y[i][j]);
			for(int j = 0; j < left_n; ++j)
				Y_ret[i].push_back(previous_y[i][j + left_l] + Y[i][j + left_window_start]);
			for(int j = left_window_end; j < n / 2; ++j)
				Y_ret[i].push_back(Y[i][j]);
		}

		/*
		printf("assign result\n");
		fflush(stdout);
		*/
		for(int i = 0; i < id.audio_channels; ++i) {
			previous_y[i].clear();
			previous_y[i].assign(Y[i].begin(), Y[i].end());
		}
	}
};

#endif
