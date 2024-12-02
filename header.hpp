#ifndef HEADER_HPP
#define HEADER_HPP

#include <vector>
#include <stack>
#include <stdlib.h>
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

struct codebook {
	unsigned int dimensions;
	unsigned int entries;
	vector<unsigned int> length;
	vector<unsigned int> codeword;
	bool ordered;
	unsigned int lookup_type;
	vector<unsigned int> multiplicands;
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
	}

	vector<unsigned int> lookup(io_buf &in) {
		for(int i = 0; i < dimension; ++i) {
		}
	}

	void setup_huffman() {
		unsigned int s = 0;
		unsigned int current_length;
		s <<= 1;
		for(unsigned int i = 0; i < entries; ++i){
			while(curreent_length < length[i]) {
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
		if(ordered == 0) {
			bool sparse = in.read_u(1);
			for(unsigned int i = 0; i < entries; ++i) {
				if(sparse == 1) {
					bool flag = in.read_u(1);
					if(flag == 1) length.push_back(in.read_u(5) + 1);
					else length.push_back(0);
				}
				else {
					length.push_back(in.read_u(5) + 1);
				}
			}
		}
		else {
			unsigned int current_entry = 0;
			unsigned int current_length = in.read_u(5) + 1;
			while(1) {
				unsigned int number = in.read_u(ilog(entries - current_entry));
				for(unsigned int i = 0; i < number; ++i)
					length.push_back(current_length);
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

		for(unsigned int i = 0; i < lookup_values; ++i)
			multiplicands.push_back(in.read_u(value_bits));

		return;
	}
};

struct floors {
	// type 0 //
	unsigned int order;
	unsigned int rate;
	unsigned int bark_map_size;
	unsigned int amplitude_bits;
	unsigned int amplitude_offset;
	unsigned int number_of_books;
	vector<int> book_list;

	void decode (io_buf &in, unsigned int type, vector<codebook> &codebook_configuration) {
		order = in.read_u(8);
		rate = in.read_u(16);
		bark_map_size = in.read_u(16);
		amplitude_bits = in.read_u(6);
		amplitude_offset = in.read_u(8);
		number_of_books = in.read_u(4) + 1;
		book_list.resize(number_of_books);
		for(unsigned int i = 0; i < number_of_books; ++i)
			book_list[i] = in.read_u(8);

		unsigned int amplitude = in.read_u(amplitude_bits);
		if(amplitude > 0) {
			vector<unsigned int> coefficients;
			unsigned int booknumber = in.read_u(ilog(number_of_books));
			if(booknumber > codebook_configuration.size()) exit(-1);
			unsigned int last = 0;
			while(1) {
				vector<unsigned int> tmp;

			}
		}
		return;
	}
};

struct setup {
	unsigned int codebook_count;
	vector<codebook> codebook_configuration;
	unsigned int time_count;
	vector<unsigned int> time;
	unsigned int floor_count;
	vector<unsigned int> floor_type;
	vector<floors> floor_configuration;

	void init(io_buf &in) {
		codebook_count = in.read_u(8) + 1;
		codebook_configuration.resize(codebook_count);
		for(unsigned int i = 0; i < codebook_count; ++i) {
			codebook_configuration[i].decode(in);
		}

		time_count = in.read_u(6) + 1;
		time.resize(time_count);
		for(unsigned int i = 0; i < time_count; ++i) {
			time[i] = in.read_u(16);
			if(time[i] != 0) exit(-1);
		}

		floor_count = in.read_u(6);
		floor_type.resize(floor_count);
		floor_configuration.resize(floor_count);
		for(unsigned int i = 0; i < floor_count; ++i) {
			floor_type[i] = in.read_u(16);
			if(floor_type[i] > 1) exit(-1);
			floor_configuration[i].decode(in, floor_type[i], codebook_configuration);
		}
	}
};

#endif
