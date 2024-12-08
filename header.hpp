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
				if(current_length == length[i] && s == codeword[i]) {
					/*
					for(int i = 0; i < entries; ++i)
						fprintf(stderr, "%lld ", codeword[i]);
					fprintf(stderr, "\n");
					*/
					//fprintf(stderr, "(entry %d code %lld length %d) ", i, s, length[i]);
					return i;
				}
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
				else length[i] = in.read_u(5) + 1;
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


#endif
