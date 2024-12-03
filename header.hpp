#ifndef HEADER_HPP
#define HEADER_HPP

#include <vector>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include "io.hpp"
#include "util.hpp"

using namespace std;

struct identification {
	unsigned int  verbis_version;
	unsigned char audio_channels;
	unsigned int  audio_sample_rate;
	int			  bitrate_maximum;
	int			  bitrate_nominal;
	int			  bitrate_minimum;
	unsigned char blocksize_0;
	unsigned char blocksize_1;
	bool		  framing_flag;

	void init (io_buf &in) {
		verbis_version    = in.read_u(32);
		audio_channels    = in.read_u(8);
		audio_sample_rate = in.read_u(32); 
		bitrate_maximum   = in.read_s(32);
		bitrate_nominal   = in.read_s(32);
		bitrate_minimum   = in.read_s(32);
		blocksize_0       = in.read_u(4);
		blocksize_1       = in.read_u(4);
		framing_flag      = in.read_u(1);
		in.read_u(7);
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
		in.read_u(7);
	}

};

struct codebooks {
	unsigned int dimensions;
	unsigned int entries;
	unsigned int* length;
	unsigned int* codeword;
	bool ordered;
	unsigned int lookup_type;
	unsigned int* multiplicands;
	float minimum_value;
	float delta_value;
	unsigned int value_bits;
	bool sequence_p;
	unsigned int lookup_values;

	unsigned int find(io_buf &in) {
		int ok = 0;
		unsigned int s = 0;
		unsigned int current_length = 0;
		while(!ok) {
			s <<= 1, s += in.read_u(1), current_length++;
			for(unsigned int i = 0; i < entries; ++i) { 
				if(current_length == length[i] && s == codeword[i]) return i;
			}
		}
		return 0;
	}

	vector<float> lookup(io_buf &in) {
		unsigned int lookup_offset = find(in);
		vector<float> value_vector(dimensions, 0);
		if(lookup_type == 1) {
			float last = 0;
			unsigned int index_divisor = 1;
			for(unsigned int i = 0; i < dimensions; ++i) {
				unsigned int multiplicand_offset = (lookup_offset / index_divisor) % lookup_values;
				value_vector[i] = multiplicands[multiplicand_offset] * delta_value + minimum_value + last;
				if(sequence_p == 1) last = value_vector[i];
				index_divisor *= lookup_values;
			}
		}
		if(lookup_type == 2) {
			float last = 0;
			unsigned int multiplicand_offset = lookup_offset * dimensions;
			for(unsigned int i = 0; i < dimensions; ++i) {
				value_vector[i] = multiplicands[multiplicand_offset] * delta_value + minimum_value + last;
				if(sequence_p == 1) last = value_vector[i];
				multiplicand_offset++;
			}
		}
		return value_vector;
	}

	vector<float> lookup(io_buf &in, unsigned int lookup_offset) {
		vector<float> value_vector(dimensions, 0);
		if(lookup_type == 1) {
			float last = 0;
			unsigned int index_divisor = 1;
			for(unsigned int i = 0; i < dimensions; ++i) {
				unsigned int multiplicand_offset = (lookup_offset / index_divisor) % lookup_values;
				value_vector[i] = multiplicands[multiplicand_offset] * delta_value + minimum_value + last;
				if(sequence_p == 1) last = value_vector[i];
				index_divisor *= lookup_values;
			}
		}
		if(lookup_type == 2) {
			float last = 0;
			unsigned int multiplicand_offset = lookup_offset * dimensions;
			for(unsigned int i = 0; i < dimensions; ++i) {
				value_vector[i] = multiplicands[multiplicand_offset] * delta_value + minimum_value + last;
				if(sequence_p == 1) last = value_vector[i];
				multiplicand_offset++;
			}
		}
		return value_vector;
	}

	void setup_huffman() {
		unsigned int s = 0;
		unsigned int current_length = 1;
		s <<= 1;
		for(unsigned int i = 0; i < entries; ++i){
			while(current_length < length[i]) {
				s <<= 1, current_length++;
			}
			while(current_length > length[i]) {
				s >>= 1, current_length--;
			}
			codeword[i] = s;
			s++;
		}
	}

	void decode(io_buf &in) {
		unsigned int tmp = in.read_u(24);
		if(tmp != 0x564342) exit(-1);
		dimensions = in.read_u(16);
		entries    = in.read_u(24);
		ordered    = in.read_u(1);
		length = new unsigned int[entries];
		codeword = new unsigned int[entries];
		if(ordered == 0) {
			bool sparse = in.read_u(1);
			for(unsigned int i = 0; i < entries; ++i) {
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
			unsigned int current_entry = 0;
			unsigned int current_length = in.read_u(5) + 1;
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
		setup_huffman();

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

		multiplicands = new unsigned int[lookup_values];
		for(unsigned int i = 0; i < lookup_values; ++i)
			multiplicands[i] = in.read_u(value_bits);

		return;
	}
};

struct floors {
	unsigned int type;
	// type 0 //
	unsigned int order;
	unsigned int rate;
	unsigned int bark_map_size;
	unsigned int amplitude_bits;
	unsigned int amplitude_offset;
	unsigned int number_of_books;
	unsigned int* book_list;

	// type 1 //
	unsigned int partitions;
	int maximum_class;
	unsigned int* partition_class_list;
	unsigned int* class_dimensions;
	unsigned int* class_subclasses;
	unsigned int* class_masterbooks;
	int** subclass_books;
	unsigned int multiplier;
	unsigned int rangebits;
	int* X_list;
	int values;

	void header_decode (io_buf &in, unsigned int _type) {
		type = _type;
		if(type == 0) {
			order = in.read_u(8);
			rate = in.read_u(16);
			bark_map_size = in.read_u(16);
			amplitude_bits = in.read_u(6);
			amplitude_offset = in.read_u(8);
			number_of_books = in.read_u(4) + 1;
			book_list = new unsigned int[number_of_books];
			for(unsigned int i = 0; i < number_of_books; ++i)
				book_list[i] = in.read_u(8);
		}
		else if(type == 1) {
			partitions = in.read_u(5);
			maximum_class = -1;
			partition_class_list = new unsigned int[partitions];
			for(unsigned int i = 0; i < partitions; ++i) {
				unsigned int tmp = in.read_u(4);
				if(maximum_class < (int) tmp) maximum_class = tmp;
				partition_class_list[i] = tmp;
			}
			maximum_class++;
			class_dimensions = new unsigned int[maximum_class];
			class_subclasses = new unsigned int[maximum_class];
			class_masterbooks = new unsigned int[maximum_class];
			subclass_books = new int*[maximum_class];
			for(int i = 0; i < maximum_class; ++i) {
				class_dimensions[i] = in.read_u(3) + 1;
				class_subclasses[i] = in.read_u(2);
				if(class_subclasses[i] != 0) class_masterbooks[i] = in.read_u(8);
				unsigned int p2 = pow(2, class_subclasses[i]);
				subclass_books[i] = new int[p2];
				for(unsigned int j = 0; j < p2; ++j) {
					subclass_books[i][j] = ((int) in.read_u(8)) - 1;
				}
			}
			multiplier = in.read_u(2) + 1;
			rangebits = in.read_u(4);
			unsigned int total = 0;
			for(unsigned int i = 0; i < partitions; ++i)
				total += class_dimensions[partition_class_list[i]];
			X_list = new int[total + 2];
			X_list[0] = 0;
			X_list[1] = pow(2, rangebits);
			values = 2;
			for(unsigned int i = 0; i < partitions; ++i) {
				unsigned int current_class_number = partition_class_list[i];
				for(unsigned int j = 0; j < class_dimensions[current_class_number]; ++j) {
					X_list[values] = in.read_u(rangebits);
					values++;
				}
			}
		}
		return;
	}
};

struct residues {
	unsigned int type;
	unsigned int begin;
	unsigned int end;
	unsigned int partition_size;
	unsigned int classification;
	unsigned int classbook;
	unsigned int* cascade;
	int** books;

	void header_decode(io_buf &in, int _type) {
		type = _type;
		begin = in.read_u(24);
		end = in.read_u(24);
		partition_size = in.read_u(24) + 1;
		classification = in.read_u(6) + 1;
		classbook = in.read_u(8);
		cascade = new unsigned int[classification];
		for(unsigned int i = 0; i < classification; ++i) {
			unsigned int high_bits = 0;
			unsigned int low_bits = in.read_u(3);
			bool bitflag = in.read_u(1);
			if(bitflag == 1) high_bits = in.read_u(5);
			cascade[i] = high_bits * 8 + low_bits;
		}
		books = new int*[classification];
		for(unsigned int i = 0; i < classification; ++i) {
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
	unsigned int coupling_step;
	unsigned int submap;
	int* magnitude;
	int* angle;
	unsigned int* mux;
	unsigned int* submap_floor;
	unsigned int* submap_residue;

	void header_decode (io_buf &in, identification &id, unsigned int fl, unsigned int re) {
		flag = in.read_u(1);
		if(flag == 1) submap = in.read_u(4) + 1;
		else submap = 1;
		flag = in.read_u(1);
		if(flag == 1) {
			coupling_step = in.read_u(8) + 1;
			magnitude = new int[coupling_step];
			angle = new int[coupling_step];
			for(unsigned int i = 0; i < coupling_step; ++i) {
				magnitude[i] = in.read_u(ilog(id.audio_channels - 1));
				angle[i] = in.read_u(ilog(id.audio_channels - 1));
				if(angle[i] == magnitude[i] 
						|| magnitude[i] > id.audio_channels - 1 
						|| angle[i] > id.audio_channels - 1) exit(-1);
			}
		}
		else coupling_step = 0;
		unsigned int rev = in.read_u(2);
		if(rev != 0) exit(-1);
		if(submap > 1) {
			mux = new unsigned int[id.audio_channels];
			for(unsigned int i = 0; i < id.audio_channels; ++i) {
				mux[i] = in.read_u(4);
				if(mux[i] > submap - 1) exit(-1);
			}
		}
		else {
			mux = new unsigned int[id.audio_channels];
			for(unsigned int i = 0; i < id.audio_channels; ++i) {
				mux[0] = 0;
			}
		}
		submap_floor = new unsigned int[submap];
		submap_residue = new unsigned int [submap];
		for(unsigned int i = 0; i < submap; ++i) {
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
	unsigned int windowtype;
	unsigned int transformtype;
	unsigned int mapping;

	void header_decode (io_buf &in, unsigned int ma) {
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
	unsigned int codebook_count;
	codebooks* codebook_configuration;
	unsigned int time_count;
	unsigned int* time;
	unsigned int floor_count;
	unsigned int* floor_type;
	floors* floor_configuration;
	unsigned int residue_count;
	unsigned int* residue_type;
	residues* residue_configuration;
	unsigned int mapping_count;
	mappings* mapping_configuration;
	unsigned int mode_count;
	modes* mode_configuration;

	void init(io_buf &in, identification &id) {
		codebook_count = in.read_u(8) + 1;
		codebook_configuration = new codebooks[codebook_count];
		for(unsigned int i = 0; i < codebook_count; ++i) {
			codebook_configuration[i].decode(in);
		}

		time_count = in.read_u(6) + 1;
		time = new unsigned int[time_count];
		for(unsigned int i = 0; i < time_count; ++i) {
			time[i] = in.read_u(16);
			if(time[i] != 0) exit(-1);
		}

		floor_count = in.read_u(6) + 1;
		floor_type = new unsigned int[floor_count];
		floor_configuration = new floors[floor_count];
		for(unsigned int i = 0; i < floor_count; ++i) {
			floor_type[i] = in.read_u(16);
			if(floor_type[i] > 1) exit(-1);
			floor_configuration[i].header_decode(in, floor_type[i]);
		}

		residue_count = in.read_u(6) + 1;
		residue_type = new unsigned int[residue_count];
		residue_configuration = new residues[residue_count];
		for(unsigned int i = 0; i < residue_count; ++i) {
			residue_type[i] = in.read_u(16);
			if(residue_type[i] > 2) exit(-1);
			residue_configuration[i].header_decode(in, residue_type[i]);
		}

		mapping_count = in.read_u(6) + 1;
		mapping_configuration = new mappings[mapping_count];
		for(unsigned int i = 0; i < mapping_count; ++i) {
			unsigned int mapping_type = in.read_u(16);
			if(mapping_type != 0) exit(-1);
			mapping_configuration[i].header_decode(in, id, floor_count, residue_count);
		}

		mode_count = in.read_u(6) + 1;
		mode_configuration = new modes[mode_count];
		for(unsigned int i = 0 ; i < mode_count; ++i) {
			mode_configuration[i].header_decode(in, mapping_count);
		}
	}

};

struct packet {
	//unsigned int type;
	unsigned int mode_number;
	unsigned int mapping_number;
	unsigned int n;
	unsigned int previous_window_flag;
	unsigned int next_window_flag;
	float *window;
	unsigned int window_center;
	unsigned int left_window_start;
	unsigned int left_window_end;
	unsigned int left_n;
	unsigned int right_window_start;
	unsigned int right_window_end;
	unsigned int right_n;

	void decode(io_buf &in, identification &id, setup &s) {
		//packet type, mode and window decode
		mode_number = in.read_u(ilog(s.mode_count - 1));
		auto &mode = s.mode_configuration[mode_number];
		mapping_number = mode.mapping;
		auto &mapping = s.mapping_configuration[mapping_number];

		unsigned int blockflag = mode.blockflag;
		n = blockflag == 0 ? id.blocksize_0 : id.blocksize_1;
		if(blockflag == 1) {
			previous_window_flag = in.read_u(1);
			next_window_flag = in.read_u(1);
		}

		window = new float[n];
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
		for(unsigned int i = 0; i < left_window_start; ++i) window[i] = 0;
		for(unsigned int i = left_window_start; i < left_window_end; ++i)
			window[i] = sin(M_PI / 2 * pow(sin((((float)i) - ((float)left_window_start) + 0.5) / left_n * ( M_PI / 2)), 2));
		for(unsigned int i = right_window_start; i < right_window_end; ++i)
			window[i] = sin(M_PI / 2 * pow(sin((((float)i) - ((float)right_window_start) + 0.5) / right_n * ( M_PI / 2) + M_PI / 2), 2));
		for(unsigned int i = right_window_end; i < n; ++i) window[i] = 0;

		// floor curve decode
		vector<vector<double>> floor_output;
		floor_output.resize(id.audio_channels);
		for(unsigned int t = 0; t < id.audio_channels; ++t)
			floor_output[t].resize(n);

		bool no_residue[id.audio_channels] = {};
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
			if(fl.type == 0) {

				unsigned int amplitude = in.read_u(fl.amplitude_bits);
				if(amplitude > 0) {
					unsigned int booknumber = in.read_u(ilog(fl.number_of_books));
					if(booknumber > s.codebook_count) exit(-1);
					float last = 0;
					vector<float> coefficients;
					while(1) {
						vector<float> tmp = s.codebook_configuration[fl.book_list[booknumber]].lookup(in);
						for(unsigned int i = 0; i < tmp.size(); ++i) tmp[i] += last;
						last = tmp[tmp.size() - 1];
						coefficients.insert(coefficients.end(), tmp.begin(), tmp.end());
						if(coefficients.size() >= fl.order) break;
					}
					/*
					// curve computation
					auto foobar = [&] (int _i) -> int {
						int ret = (int) (bark(((double) fl.rate) * _i / 2 / n) * fl.bark_map_size / bark(0.5 * fl.rate));
						return ret;
					};
					vector<float> map(n, 0);
					for(unsigned int i = 0; i < n - 1; ++i) {
						map[i] = min((int) (fl.bark_map_size - 1), foobar(i));
					}
					map[n - 1] = -1;

					auto getpq = [&] (double omega, double &p, double &q) -> void {
						double c = cos(omega);
						p = 1, q = 1;
						if((fl.order & 1) == 0) {
							for(unsigned int j = 0; j <= ((fl.order - 3) / 2); ++j)
								p *= (coefficients[2 * j + 1] - c) * (coefficients[2 * j + 1 - c] - c);
							for(unsigned int j = 0; j <= ((fl.order - 1) / 2); ++j)
								q *= (coefficients[2 * j] - c) * (coefficients[2 * j] - c);
							p *= (1 - c * c) * 4;
							q /= 4;
						}
						else {
							for(unsigned int j = 0; j <= ((fl.order - 2) / 2); ++j)
								p *= (coefficients[2 * j + 1] - c) * (coefficients[2 * j + 1 - c] - c);
							for(unsigned int j = 0; j <= ((fl.order - 2) / 2); ++j)
								q *= (coefficients[2 * j] - c) * (coefficients[2 * j] - c);
							p *= (1 - c) * 2;
							q *= (1 + c) * 2;
						}
					};

					auto linear = [&] (double p, double q) -> double {
						return exp(.11512925 * (((double) amplitude) * fl.amplitude_offset / (pow(2, fl.amplitude_bits) - 1) / sqrt(p + q) - fl.amplitude_offset));
					};

					for(unsigned int i = 0; i < n; ) {
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
					*/
				}
				else {
					unused = 1;
					for(unsigned int i = 0; i < n; ++i)
						floor_output[t][i] = 0;
				}
			}

			else {
				bool nonzero = in.read_u(1);
				if(nonzero == 0) {
					unused = 1;
				}
				else {
					unsigned int rr[4] = {256, 128, 86, 64};
					unsigned int range = rr[fl.multiplier - 1];
					unsigned int Y[2 + fl.values];
					Y[0] = in.read_u(ilog(range - 1));
					Y[1] = in.read_u(ilog(range - 1));
					unsigned int offset = 2;
					for(unsigned int i = 0; i < fl.partitions; ++i) {
						unsigned int class_ = fl.partition_class_list[i];
						unsigned int cdim = fl.class_dimensions[class_];
						unsigned int cbits = fl.class_subclasses[class_];
						unsigned int csub = pow(2, cbits) - 1;
						unsigned int cval = 0;
						if(cbits > 0) cval = s.codebook_configuration[fl.class_masterbooks[class_]].find(in);
						for(unsigned int j = 0; j < cdim; ++j) {
							unsigned int book = fl.subclass_books[class_][cval & csub];
							cval >>= cbits;
							if(book >= 0) Y[j + offset] = s.codebook_configuration[book].find(in);
							else Y[j + offset] = 0;
						}
						offset += cdim;
					}
					/*
					// step 1
					int final_Y[range];
					bool step2_flag[2 + fl.values];
					step2_flag[0] = 1, step2_flag[1] = 1;
					final_Y[0] = Y[0], final_Y[1] = Y[1];
					for(int i = 2; i < fl.values; ++i) {
						int low_neighbor_offset = low_neighbor(fl.X_list, i);
						int high_neighbor_offset = high_neighbor(fl.X_list, i, fl.values + 2);
						int predicted = render_point(
								fl.X_list[low_neighbor_offset], final_Y[low_neighbor_offset],
								fl.X_list[high_neighbor_offset], final_Y[high_neighbor_offset], fl.X_list[i]);
						int val = Y[i];
						int highroom = range - predicted, lowroom = predicted, room;
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
					floor1_sort(fl.X_list, final_Y, step2_flag, 2 + fl.values);
					int hx = 0, lx = 0;
					int ly = final_Y[0] * fl.multiplier, hy = 0;
					int ma_x = 0;
					for(int i = 0; i < fl.values; ++i) ma_x = max(ma_x, fl.X_list[i]);
					vector<int> output(ma_x, 0);
					for(int i = 1; i < fl.values; ++i) {
						if(step2_flag[i] == 1) {
							hy = final_Y[i] * fl.multiplier;
							hx = fl.X_list[i];
							render_line(lx, ly, hx, hy, output);
							lx = hx, ly = hy;
						}
					}
					if(hx < (int) n) render_line(hx, hy, n, hy, output);
					if(hx > (int) n) output.resize(n);
					for(unsigned int i = 0; i < n; ++i) {
						floor_output[t][i] = floor1_inverse_dB_table[output[i]];
					}
					*/
				}
			}
			if(unused) {
				for(unsigned int i = 0; i < n; ++i) floor_output[t][i] = 0;
				no_residue[t] = 1;
			}
		}

		// nonzero vector propagate
		vector<int> residue_output(mapping.submap);
		for(unsigned int i = 0; i < mapping.submap; ++i) {
			unsigned int ch = 0;
			bool do_not_decode_flag[id.audio_channels] = {};
			for(unsigned int j = 0; j < id.audio_channels; ++j) {
				if(mapping.mux[j] == i) {
					if(no_residue[j] == 1) do_not_decode_flag[ch] = 1;
					else do_not_decode_flag[ch] = 0;
					ch++;
				}
			}
			unsigned int residue_number = mapping.submap_residue[i];
			auto &residue = s.residue_configuration[residue_number];
			unsigned int residue_type = s.residue_type[residue_number];
			// decode packet
			unsigned int actual_size = n / 2;
			if(residue_type == 2) actual_size *= ch;
			unsigned int limit_residue_begin = min(residue.begin, actual_size);
			unsigned int limit_residue_end = min(residue.end, actual_size);
			unsigned int classwords_per_codeword = s.codebook_configuration[residue.classbook].dimensions;
			unsigned int n_to_read = limit_residue_end - limit_residue_begin;
			unsigned int partitions_to_read = n_to_read / residue.partition_size;
			vector<vector<float>> r_output(ch);

			if(n_to_read != 0) {
				for(int pass = 0; pass < 8; ++pass) {
					unsigned int partition_count = 0;
					while(partition_count < partitions_to_read) {
						int classifications[ch][classwords_per_codeword + partition_count];
						if(pass == 0) {
							for(unsigned int j = 0; j < ch; ++j) if(do_not_decode_flag[j] == 0) {
								unsigned int temp = s.codebook_configuration[residue.classbook].find(in);
								for(int t = classwords_per_codeword - 1; t >= 0; --t) {
									classifications[j][t + partition_count] = temp % residue.classification;
									temp /= residue.classification;
								}
							}
						}
						for(int t = 0; t < classwords_per_codeword && partition_count < partitions_to_read; ++t) {
							for(unsigned int j = 0; j < ch; ++j) if(do_not_decode_flag[j] == 0) {
								int vqclass = classifications[j][partition_count];
								int vqbook = residue.books[vqclass][pass];
								if(vqbook != -1) {
									r_output[t] = s.codebook_configuration[vqbook].lookup(in, limit_residue_begin + partition_count * residue.partition_size);
								}
							}
							partition_count++;
						}
					}
				}
			}
		}
	}
};

#endif
