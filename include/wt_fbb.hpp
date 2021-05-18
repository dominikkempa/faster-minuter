/**
 * @file    wt_fbb.cpp
 * @section LICENCE
 *
 * This file is part of the "Faster Minuter" index v0.1.0
 * See: https://github.com/dominikkempa/faster-minuter
 *
 * Copyright (C) 2015-2021
 *   Juha Karkkainen <juha.karkkainen (at) cs.helsinki.fi>
 *   Dominik Kempa <dominik.kempa (at) gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 **/

#ifndef __WT_FBB_HPP_INCLUDED
#define __WT_FBB_HPP_INCLUDED

#include <cstdint>
#include <vector>
#include <queue>
#include <algorithm>

#include <sdsl/sdsl_concepts.hpp>
#include <sdsl/int_vector_buffer.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/io.hpp>


#define ADD_NAVIGATIONAL_BLOCK_HEADER  // minimally increases space but removes navigational queries, default: ON
//#define SPARSE_SUPERBLOCK_MAPPING    // reduces the space at the cost of rank query slowdown, default: OFF
#define ALLOW_VARIABLE_BLOCK_SIZE      // permits different blocks sizes, reduces time and space, default: ON
#define FAST_CONSTRUCTION              // use fast computation of block sizes, default: ON


//! Wavelet tree based on a fixed block boosting technique.
/*!
 *  \tparam t_bitvector   Underlying bitvector structure.
 *  \tparam t_rank        Rank support for pattern `1` on the bitvector.
 *  \tparam t_sbs_log     Base 2 logarithm of the superblock size.
 *
 *  @ingroup wt
 *
 *  \par References:
 *     [1] Simon Gog, Juha Karkkainen, Dominik Kempa,
 *         Matthias Petri, Simon J. Puglisi:
 *         Faster, Minuter.
 *         DCC 2016: 53-62
 *     [2] Simon Gog, Juha Kärkkäinen, Dominik Kempa,
 *         Matthias Petri, Simon J. Puglisi:
 *         Fixed Block Compression Boosting in FM-Indexes: Theory
 *         and Practice.
 *         Algorithmica 81(4): 1370-1391 (2019)
 */
template<class t_bitvector       = sdsl::hyb_vector<>,
         class t_rank            = typename t_bitvector::rank_1_type,
#ifndef ALLOW_VARIABLE_BLOCK_SIZE
         std::uint64_t t_bs_log  = 14,
#endif
         std::uint64_t t_sbs_log = 20>
class wt_fbb {
  public:
    typedef sdsl::wt_tag index_category;
    typedef std::uint64_t size_type;
    using alphabet_category = sdsl::byte_alphabet_tag;

  private:
#ifndef ALLOW_VARIABLE_BLOCK_SIZE
    static const std::uint64_t k_block_size_log;
    static const std::uint64_t k_block_size;
#endif
    static const std::uint64_t k_superblock_size_log;
    static const std::uint64_t k_superblock_size;
    static const std::uint64_t k_hyperblock_size;

    struct block_header_item {
#ifdef ADD_NAVIGATIONAL_BLOCK_HEADER
      std::uint32_t m_bv_rank;                 // rank at the beginning of block vector
#endif
      std::uint32_t m_bv_offset;               // offset in the superblock bitvector
      std::uint32_t m_var_size_header_offset;  // offset in the array storing variable-size block header
      std::uint8_t  m_sigma;                   // block alphabet size - 1
      std::uint8_t  m_tree_height;             // height of the wavelet tree for block
    } __attribute__ ((packed));

    struct superblock_header_item {
      std::uint8_t m_sigma;                            // superblock alphabet size - 1
#ifdef ALLOW_VARIABLE_BLOCK_SIZE
      std::uint8_t m_block_size_log;
#endif
      t_bitvector m_bitvector;                         // superblock bitvector
      t_rank m_rank_support;                           // rank support for superblock bitvector
      std::vector<std::uint8_t> m_var_block_headers;   // variable-size block header
      std::vector<block_header_item> m_block_headers;  // fixed-size block header
#ifdef SPARSE_SUPERBLOCK_MAPPING
      std::vector<std::uint8_t> m_mapping_data;
      sdsl::bit_vector_il<> m_mapping_bv;
      typename sdsl::bit_vector_il<>::rank_1_type m_mapping_bv_rank_support;
      typename sdsl::bit_vector_il<>::select_1_type m_mapping_bv_select_support;
      std::uint32_t m_one_bit_count;
#else
      std::vector<std::uint8_t> m_mapping;            // mapping from superblock alphabet to block alphabet
#endif
    };

    friend std::uint64_t serialize(const superblock_header_item &shi, std::ostream &out, sdsl::structure_tree_node* v = nullptr, std::string name = "") {
      sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(shi));
      std::uint64_t written_bytes = 0;
      written_bytes += sdsl::serialize(shi.m_sigma, out, child, "sigma");
#ifdef ALLOW_VARIABLE_BLOCK_SIZE
      written_bytes += sdsl::serialize(shi.m_block_size_log, out, child, "block_size_log");
#endif
      written_bytes += sdsl::serialize(shi.m_bitvector, out, child, "bitvector");
      written_bytes += sdsl::serialize(shi.m_rank_support, out, child, "rank_support");
      written_bytes += sdsl::serialize(shi.m_var_block_headers, out, child, "var_block_headers");
      written_bytes += sdsl::serialize(shi.m_block_headers, out, child, "block_headers");
#ifdef SPARSE_SUPERBLOCK_MAPPING
      written_bytes += sdsl::serialize(shi.m_mapping_data, out, child, "mapping_data");
      written_bytes += sdsl::serialize(shi.m_mapping_bv, out, child, "mapping_bv");
      written_bytes += sdsl::serialize(shi.m_mapping_bv_rank_support, out, child, "rank_support");
      written_bytes += sdsl::serialize(shi.m_mapping_bv_select_support, out, child, "select_support");
      written_bytes += sdsl::serialize(shi.m_one_bit_count, out, child, "one_bit_count");
#else
      written_bytes += sdsl::serialize(shi.m_mapping, out, child, "mapping");
#endif
      sdsl::structure_tree::add_size(child, written_bytes);
      return written_bytes;
    }

    friend void load(superblock_header_item &shi, std::istream &in) {
      sdsl::load(shi.m_sigma, in);
#ifdef ALLOW_VARIABLE_BLOCK_SIZE
      sdsl::load(shi.m_block_size_log, in);
#endif
      sdsl::load(shi.m_bitvector, in);
      sdsl::load(shi.m_rank_support, in);
      sdsl::load(shi.m_var_block_headers, in);
      sdsl::load(shi.m_block_headers, in);
#ifdef SPARSE_SUPERBLOCK_MAPPING
      sdsl::load(shi.m_mapping_data);
      sdsl::load(shi.m_mapping_bv);
      sdsl::load(shi.m_rank_support);
      sdsl::load(shi.m_select_support);
      sdsl::load(shi.m_one_bit_count);
      shi.m_rank_support.set_vector(&shi.m_mapping_bv);
      shi.m_select_support.set_vector(&shi.m_mapping_bv);
#else
      sdsl::load(shi.m_mapping, in);
#endif
      shi.m_rank_support.set_vector(&shi.m_bitvector);
    }

  private:
    std::uint64_t m_size;                                      // original sequence length
    std::vector<std::uint64_t> m_count;                        // global symbol counts
    std::vector<std::uint64_t> m_hyperblock_rank;              // ranks at hyperblock boundary
    std::vector<std::uint32_t> m_superblock_rank;              // ranks at superblock boundary
    std::vector<std::uint8_t> m_global_mapping;                // mapping from global alphabet to superblock alphabet
    std::vector<superblock_header_item> m_superblock_headers;  // superblock headers

    void copy(const wt_fbb& tree) {
      m_size = tree.m_size;
      m_count = tree.m_count;
      m_hyperblock_rank = tree.m_hyperblock_rank;
      m_superblock_rank = tree.m_superblock_rank;
      m_global_mapping = tree.m_global_mapping;
      m_superblock_headers = tree.m_superblock_headers;
    }

  public:
    // Default constructor
    wt_fbb() = default;

    // Copy constructor
    wt_fbb(const wt_fbb& tree) {
      copy(tree);
    }

    // Move constructor
    wt_fbb(wt_fbb&& tree) {
      *this = std::move(tree);
    }

    // Constructor
    wt_fbb(sdsl::int_vector_buffer<(std::uint8_t)8> &text_buf, std::uint64_t text_length) {
      init(text_buf, text_length);
    }

    // Constructor
    wt_fbb(const std::uint8_t *text, std::uint64_t text_length) {
      init(text, text_length);
    }

  private:
    static void compute_symbol_freq(
        const std::uint8_t *text,
        std::uint64_t text_length,
        std::vector<std::uint64_t> &freq) {

      std::fill(freq.begin(), freq.end(), 0UL);
      for (std::uint64_t i = 0; i < text_length; ++i)
        ++freq[text[i]];
    }

    static void compute_huffman_code_lengths(
        std::vector<std::uint64_t> &freq,
        std::vector<std::uint64_t> &codelen) {

      std::fill(codelen.begin(), codelen.end(), 0UL);
      typedef std::pair<std::uint64_t, std::vector<std::uint8_t> > pq_item;
      std::priority_queue<pq_item,
        std::vector<pq_item>, std::greater<pq_item> > pq;

      for (std::uint64_t i = 0; i < 256; ++i) {
        if (freq[i] > 0) {
          std::vector<std::uint8_t> v { (std::uint8_t)i };
          pq.push(std::make_pair(freq[i], v));
        }
      }

      while (pq.size() > 1) {
        pq_item x = pq.top(); pq.pop();
        pq_item y = pq.top(); pq.pop();
        std::vector<std::uint8_t> v = x.second;
        v.insert(v.end(), y.second.begin(), y.second.end());
        for (std::uint8_t i : v) ++codelen[i];
        pq.push(std::make_pair(x.first + y.first, v));
      }
    }

    static void assign_canonical_huffman_codes(
        std::vector<std::uint64_t> &freq,
        std::vector<std::uint64_t> &codelen,
        std::vector<std::uint64_t> &code) {

      std::fill(code.begin(), code.end(), 0UL);
      std::vector<std::pair<std::uint64_t, std::uint8_t> > sym;
      for (std::uint64_t i = 0; i < 256; ++i)
        if (freq[i] > 0)
          sym.push_back(std::make_pair(codelen[i], (std::uint8_t)i));
      std::sort(sym.begin(), sym.end());

      for (std::uint64_t c = 0, i = 0; i < sym.size(); ++i) {
        if (i == 0) code[sym[i].second] = c;
        else {
          c = (c + 1) << (sym[i].first - sym[i - 1].first);
          code[sym[i].second] = c;
        }
      }
    }

    // Restore the Huffman code of the given symbol in block alphabet
    // (`block_c') from the block header and return via variables `code'
    // and `codelen'. We assume tree_height > 0.
    static void restore_code_from_block_header(
        std::uint64_t block_c,
        const std::uint8_t *block_header_ptr,
        std::uint64_t tree_height,
        std::uint64_t &code,
        std::uint64_t &codelen) {

      code = 0;
      codelen = 1;
      std::uint64_t leaf_count = 0;
      while (codelen < tree_height) {
        code <<= 1;
        std::uint64_t this_level_leaf_count = *block_header_ptr;
        if (leaf_count + this_level_leaf_count > block_c) {
          code += (block_c - leaf_count);
          break;
        } else {
          code += this_level_leaf_count;
          ++codelen;
          leaf_count += this_level_leaf_count;
          block_header_ptr += 3;
        }
      }
      if (codelen == tree_height) {
        code <<= 1;
        code += (block_c - leaf_count);
      }
    }

    static std::uint64_t compute_symbol_from_block_header(
        const std::uint8_t *block_header_ptr,
        std::uint64_t code,
        std::uint64_t codelen) {
      std::uint64_t block_c = 0;
      std::uint64_t temp_code = 0;
      for (std::uint64_t i = 1; i < codelen; ++i) {
        std::uint64_t this_level_leaf_count = *block_header_ptr;
        block_header_ptr += 3;
        temp_code += this_level_leaf_count;
        block_c += this_level_leaf_count;
        temp_code <<= 1;
      }
      block_c += code - temp_code;
      return block_c;
    }

    void encode_block(
        const std::uint8_t *block_ptr,
        std::vector<std::uint64_t> &block_rank,
        std::uint64_t block_size,
        sdsl::bit_vector &superblock_bv,
        std::uint64_t &ones_count,
        std::uint64_t superblock_bv_offset,
        std::uint8_t *block_header_ptr) {

      // Compute Huffman code.
      std::vector<std::uint64_t> freq(256), codelen(256), code(256);
      compute_symbol_freq(block_ptr, block_size, freq);
      compute_huffman_code_lengths(freq, codelen);
      assign_canonical_huffman_codes(freq, codelen, code);
      std::uint64_t max_code_length =
        *std::max_element(codelen.begin(), codelen.end());

      ones_count = 0;
#ifdef ADD_NAVIGATIONAL_BLOCK_HEADER
      std::vector<std::uint64_t> ones_in_bv(256, 0UL);
#endif

      // Compute bitvectors for all internal nodes
      // of the tree and append to superblock_bv.
      if (std::count(freq.begin(), freq.end(), 0UL) < 255) {

        // Collect IDs of all internal nodes in the tree. The ID of the
        // node is a number whose bits are taken from the root-to-node path
        // (first bit on the path is MSB in ID) prepended with 1, e.g., the
        // node with path 011 has ID 1011 = (DEC)11. This encoding has a
        // number of nice properties: (1) every possible node has unique
        // ID, (2) sorted IDs correspond to nodes in BFS order (and
        // left-to-right within level), which is the order in which we
        // concatenate bitvectors, (3) ID of a sibling with ID x is (x xor 1).
        std::vector<std::uint64_t> internal_node_ids;
        for (std::uint64_t i = 0; i < 256; ++i) {
          if (freq[i] > 0) {
            for (std::uint64_t depth = 0; depth < codelen[i]; ++depth) {
              std::uint64_t id =
                (((1UL << codelen[i]) | code[i]) >> (codelen[i] - depth));
              internal_node_ids.push_back(id);
            }
          }
        }
        std::sort(internal_node_ids.begin(), internal_node_ids.end());
        internal_node_ids.erase(std::unique(internal_node_ids.begin(),
            internal_node_ids.end()), internal_node_ids.end());

        // Compute the mapping from internal nodes to bitvectors (which are
        // numbered with consecutive numbers starting from 0, according to
        // the order in which they are concatenated).
        std::vector<std::uint64_t>
          internal_node_bv_id(1UL << max_code_length);
        for (std::uint64_t i = 0; i < internal_node_ids.size(); ++i)
          internal_node_bv_id[internal_node_ids[i]] = i;

        // Compute the size of bitvector for every internal node.
        std::vector<std::uint64_t>
          internal_node_bv_size(internal_node_ids.size(), 0UL);
        for (std::uint64_t i = 0; i < 256; ++i) {
          if (freq[i] > 0) {
            for (std::uint64_t depth = 0; depth < codelen[i]; ++depth) {
              std::uint64_t id =
                (((1UL << codelen[i]) | code[i]) >> (codelen[i] - depth));
              internal_node_bv_size[internal_node_bv_id[id]] += freq[i];
            }
          }
        }

        // Allocate bitvectors for all internal nodes.
        std::vector<std::uint8_t> **internal_node_bv =
          new std::vector<std::uint8_t>*[internal_node_ids.size()];
        for (std::uint64_t i = 0; i < internal_node_ids.size(); ++i)
          internal_node_bv[i] =
            new std::vector<std::uint8_t>(internal_node_bv_size[i], 0);

        // Fill in the bitvectors for all internal nodes.
        std::vector<std::uint64_t>
          node_visit_count(1UL << (max_code_length + 1), 0UL);
        for (std::uint64_t i = 0; i < block_size; ++i) {
          std::uint8_t sym = block_ptr[i];
          std::uint64_t pos = i;
          for (std::uint64_t depth = 0; depth < codelen[sym]; ++depth) {
            std::uint64_t id =
              (((1UL << codelen[sym]) | code[sym]) >> (codelen[sym] - depth));
            if (depth > 0) {
              pos -= node_visit_count[id ^ 1];
              ++node_visit_count[id];
            }
            if (code[sym] & (1UL << (codelen[sym] - depth - 1))) {
              (*internal_node_bv[internal_node_bv_id[id]])[pos] = 1;
#ifdef ADD_NAVIGATIONAL_BLOCK_HEADER
              ones_in_bv[internal_node_bv_id[id]] += 1;
#endif
              ++ones_count;
            }
          }
          ++node_visit_count[(1UL << codelen[sym]) | code[sym]];
        }

        // Append bitvectors of internal nodes to superblock bitvector.
        for (std::uint64_t i = 0; i < internal_node_ids.size(); ++i)
          for (std::uint64_t j = 0; j < internal_node_bv[i]->size(); ++j)
            superblock_bv[superblock_bv_offset++] = (*internal_node_bv[i])[j];

        // Clean up.
        for (std::uint64_t i = 0; i < internal_node_ids.size(); ++i)
          delete internal_node_bv[i];
        delete[] internal_node_bv;
      }

      // Fill in the variable-size block header.
      {
        std::vector<std::uint64_t> codelen_freq(max_code_length, 0UL);
        for (std::uint64_t i = 0; i < 256; ++i)
          if (freq[i] > 0 && codelen[i] < max_code_length)
            ++codelen_freq[codelen[i]];  // encoded with 1 byte

        // level_total_freq[d] = total frequency of symbols that have code
        // length longer than d, i.e., the total length of bitvectors
        // associated with internal nodes at depth d in the tree.
        std::vector<std::uint64_t> level_total_freq(max_code_length, 0UL);
        for (std::uint64_t i = 0; i < 256; ++i)
          if (freq[i] > 0)
            for (std::uint64_t depth = 1; depth < codelen[i]; ++depth)
              level_total_freq[depth] += freq[i];  // encoded with 2 bytes

        // Store the number of leaves and total size of bitvectors
        // corresponding to internal nodes at each level in the tree
        // (except root level and deepest levels) minus one. Using
        // 2 bytes limits the block size to 2^16.
        for (std::uint64_t depth = 1; depth < max_code_length; ++depth) {
          *(block_header_ptr++) = (std::uint8_t)codelen_freq[depth];                          // leaf count, 1 byte
          *((std::uint16_t *)block_header_ptr) = (std::uint16_t)level_total_freq[depth] - 1;  // total bv size minus one, 2 bytes
          block_header_ptr += 2;
        }

        // Sort symbols by frequency.
        std::vector<std::pair<std::uint64_t, std::uint8_t> > sym;
        for (std::uint64_t i = 0; i < 256; ++i)
          if (freq[i] > 0) sym.push_back(std::make_pair(codelen[i], (std::uint8_t)i));
        std::sort(sym.begin(), sym.end());

        // Store rank value at block boundary (with respect to to superblock
        // boundary) and global symbol for each leaf. Note: using 3 bytes
        // to store rank at block boundary limits the superblock size to 2^24.
        for (std::uint64_t i = 0; i < sym.size(); ++i) {
          std::uint64_t symbol = (std::uint64_t)sym[i].second;  // encoded with 1 byte
          std::uint64_t rank_value = block_rank[symbol];        // encoded with 3 bytes
          *((std::uint32_t *)block_header_ptr) = (std::uint32_t)(symbol | (rank_value << 8));
          block_header_ptr += 4;
        }

#ifdef ADD_NAVIGATIONAL_BLOCK_HEADER
        // For every internal node, store the number of 1-bits in the
        // bitvector corresponding to that node and all its left siblings
        // (excluding leaves).
        std::uint64_t number_of_internal_nodes_current_level = 1;
        for (std::uint64_t depth = 0, ptr = 0; depth < max_code_length; ++depth) {
          std::uint64_t one_bits_current_level_count = 0;
          for (std::uint64_t j = 0; j < number_of_internal_nodes_current_level; ++j) {
            one_bits_current_level_count += ones_in_bv[ptr++];
            *((std::uint16_t *)block_header_ptr) = (std::uint16_t)one_bits_current_level_count;  // 2 bytes
            block_header_ptr += 2;
          }
          if (depth + 1 != max_code_length) {
            std::uint64_t next_level_leaf_count = codelen_freq[depth + 1];
            number_of_internal_nodes_current_level <<= 1;
            number_of_internal_nodes_current_level -= next_level_leaf_count;
          }
        }
#endif
      }
    }

    // Compute block bitvector. The result is
    // appended to superblock bitvector.
    void compute_block_bv(
        const std::uint8_t *block_ptr,
        std::uint64_t block_size,
        sdsl::bit_vector &superblock_bv,
        std::uint64_t superblock_bv_offset) {

      // Compute Huffman code.
      std::vector<std::uint64_t> freq(256), codelen(256), code(256);
      compute_symbol_freq(block_ptr, block_size, freq);
      compute_huffman_code_lengths(freq, codelen);
      assign_canonical_huffman_codes(freq, codelen, code);
      std::uint64_t max_code_length =
        *std::max_element(codelen.begin(), codelen.end());

      if (std::count(freq.begin(), freq.end(), 0UL) < 255) {

        // Collect IDs of all internal nodes in the tree.
        std::vector<std::uint64_t> internal_node_ids;
        for (std::uint64_t i = 0; i < 256; ++i) {
          if (freq[i] > 0) {
            for (std::uint64_t depth = 0; depth < codelen[i]; ++depth) {
              std::uint64_t id = 
                (((1UL << codelen[i]) | code[i]) >> (codelen[i] - depth));
              internal_node_ids.push_back(id);
            }
          }
        }
        std::sort(internal_node_ids.begin(), internal_node_ids.end());
        internal_node_ids.erase(std::unique(internal_node_ids.begin(),
            internal_node_ids.end()), internal_node_ids.end());

        // Compute the mapping from internal nodes to bitvectors (which are
        // numbered with consecutive numbers starting from 0, according to
        // the order in which they are concatenated).
        std::vector<std::uint64_t> internal_node_bv_id(1UL << max_code_length);
        for (std::uint64_t i = 0; i < internal_node_ids.size(); ++i)
          internal_node_bv_id[internal_node_ids[i]] = i;

        // Compute the size of bitvector for every internal node.
        std::vector<std::uint64_t>
          internal_node_bv_size(internal_node_ids.size(), 0UL);
        for (std::uint64_t i = 0; i < 256; ++i) {
          if (freq[i] > 0) {
            for (std::uint64_t depth = 0; depth < codelen[i]; ++depth) {
              std::uint64_t id =
                (((1UL << codelen[i]) | code[i]) >> (codelen[i] - depth));
              internal_node_bv_size[internal_node_bv_id[id]] += freq[i];
            }
          }
        }

        // Allocate bitvectors for all internal nodes.
        std::vector<std::uint8_t> **internal_node_bv =
          new std::vector<std::uint8_t>*[internal_node_ids.size()];
        for (std::uint64_t i = 0; i < internal_node_ids.size(); ++i)
          internal_node_bv[i] =
            new std::vector<std::uint8_t>(internal_node_bv_size[i], 0);

        // Fill in the bitvectors for all internal nodes.
        std::vector<std::uint64_t>
          node_visit_count(1UL << (max_code_length + 1), 0UL);
        for (std::uint64_t i = 0; i < block_size; ++i) {
          std::uint8_t sym = block_ptr[i];
          std::uint64_t pos = i;
          for (std::uint64_t depth = 0; depth < codelen[sym]; ++depth) {
            std::uint64_t id =
              (((1UL << codelen[sym]) | code[sym]) >> (codelen[sym] - depth));
            if (depth > 0) {
              pos -= node_visit_count[id ^ 1];
              ++node_visit_count[id];
            }
            if (code[sym] & (1UL << (codelen[sym] - depth - 1)))
              (*internal_node_bv[internal_node_bv_id[id]])[pos] = 1;
          }
          ++node_visit_count[(1UL << codelen[sym]) | code[sym]];
        }

        // Append bitvectors of internal nodes to superblock bitvector.
        for (std::uint64_t i = 0; i < internal_node_ids.size(); ++i)
          for (std::uint64_t j = 0; j < internal_node_bv[i]->size(); ++j)
            superblock_bv[superblock_bv_offset++] = (*internal_node_bv[i])[j];

        // Clean up.
        for (std::uint64_t i = 0; i < internal_node_ids.size(); ++i)
          delete internal_node_bv[i];
        delete[] internal_node_bv;
      }
    }

    void encode_blocks_in_superblock(
        const std::uint8_t *superblock_ptr,
        std::uint64_t superblock_id,
        std::uint64_t block_size_log) {

      superblock_header_item &superblock_header =
        m_superblock_headers[superblock_id];

      std::uint64_t block_size = (1UL << block_size_log);
      std::uint64_t superblock_beg = superblock_id * k_superblock_size;
      std::uint64_t superblock_end = std::min(superblock_beg + k_superblock_size, m_size);
      std::uint64_t this_superblock_size = superblock_end - superblock_beg;
      std::uint64_t this_superblock_sigma = (std::uint64_t)superblock_header.m_sigma + 1;
#ifdef ALLOW_VARIABLE_BLOCK_SIZE
      superblock_header.m_block_size_log = block_size_log;
#endif

#ifdef SPARSE_SUPERBLOCK_MAPPING

      // Create temporary containers for superblock mapping.
      sdsl::bit_vector this_superblock_mapping_bv(
          this_superblock_sigma * (k_superblock_size / block_size), 0);
      std::vector<std::uint8_t> this_superblock_mapping_data[256];
      std::uint32_t one_bit_count = 0;
#else

      // Allocate the array storing mapping
      // from superblock alphabet to block alphabets.
      superblock_header.m_mapping = std::vector<std::uint8_t>(
          this_superblock_sigma * (k_superblock_size / block_size), 255);
#endif

      // Fill in the fixed-size block header and superblock headers, and
      // compute sizes of variable-size block headers and bitvectors.
      std::uint64_t superblock_bv_size = 0;
      std::uint64_t variable_block_header_size = 0;
      std::uint64_t blocks_in_this_superblock = (this_superblock_size + block_size - 1) / block_size;
      superblock_header.m_block_headers = std::vector<block_header_item>(blocks_in_this_superblock);
      for (std::uint64_t block_id = 0; block_id < blocks_in_this_superblock; ++block_id) {
        block_header_item &block_header = superblock_header.m_block_headers[block_id];
        std::uint64_t block_beg = block_id * block_size;
        std::uint64_t block_end = std::min(block_beg + block_size, this_superblock_size);
        std::uint64_t this_block_size = block_end - block_beg;
        const std::uint8_t *block_ptr = superblock_ptr + block_beg;

        // Compute global-to-block symbol mapping and basic
        // information about symbol distribution in the block.
        std::uint64_t tree_height = 0, bv_size = 0, sigma = 0;
        std::vector<std::uint16_t> global_to_block_mapping(256, 256);
        {

          // Compute Huffman code lengths.
          std::vector<std::uint64_t> freq(256), codelen(256);
          compute_symbol_freq(block_ptr, this_block_size, freq);
          compute_huffman_code_lengths(freq, codelen);

          // Sort symbols by frequency.
          std::vector<std::pair<std::uint64_t, std::uint8_t> > sym;
          for (std::uint64_t i = 0; i < 256; ++i)
            if (freq[i] > 0)
              sym.push_back(std::make_pair(codelen[i], (std::uint8_t)i));
          std::sort(sym.begin(), sym.end());

          // Fill in all fields.
          sigma = sym.size();
          tree_height = *std::max_element(codelen.begin(), codelen.end());
          for (std::uint64_t i = 0; i < sym.size(); ++i) {
            global_to_block_mapping[sym[i].second] = (std::uint16_t)i;
            if (sym.size() > 1)
              bv_size += freq[sym[i].second] * codelen[sym[i].second];
          }
        }

        // Fill in the fixed-size block header.
        block_header.m_bv_offset = (std::uint32_t)superblock_bv_size;
        block_header.m_var_size_header_offset = (std::uint32_t)variable_block_header_size;
        block_header.m_tree_height = (std::uint8_t)tree_height;
        block_header.m_sigma = (std::uint8_t)(sigma - 1);

        // Store the superblock mapping.
        for (std::uint64_t i = 0; i < 256; ++i) {
          if (global_to_block_mapping[i] != 256) {
            std::uint8_t superblock_char = m_global_mapping[superblock_id * 256 + i];
            std::uint64_t addr = (std::uint64_t)superblock_char * (k_superblock_size / block_size) + block_id;
            std::uint8_t mapping_val = std::min((std::uint16_t)254, global_to_block_mapping[i]);
#ifdef SPARSE_SUPERBLOCK_MAPPING
            this_superblock_mapping_bv[addr] = 1;
            this_superblock_mapping_data[superblock_char].push_back(mapping_val);
            ++one_bit_count;
#else
            superblock_header.m_mapping[addr] = mapping_val;
#endif
          }
        }

        // Update the size of superblock bitvector.
        superblock_bv_size += bv_size;

        // Update the total size of variable-size block headers.
        if (tree_height > 1)
          variable_block_header_size += (tree_height - 1) * 3;  // info about each level in the tree
        variable_block_header_size += sigma * 4;                // info about each leaf in the tree
#ifdef ADD_NAVIGATIONAL_BLOCK_HEADER
        variable_block_header_size += (sigma - 1) * 2;          // additional navigational info
#endif
      }

#ifdef SPARSE_SUPERBLOCK_MAPPING
      superblock_header.m_mapping_bv = sdsl::bit_vector_il<>(this_superblock_mapping_bv);
      superblock_header.m_mapping_bv_rank_support = sdsl::bit_vector_il<>::rank_1_type(&superblock_header.m_mapping_bv);
      superblock_header.m_mapping_bv_select_support = sdsl::bit_vector_il<>::select_1_type(&superblock_header.m_mapping_bv);

      std::uint64_t total_sparse_mapping_size = 0;
      for (std::uint64_t i = 0; i < 256; ++i)
        total_sparse_mapping_size += this_superblock_mapping_data[i].size();
      superblock_header.m_mapping_data = std::vector<std::uint8_t>(total_sparse_mapping_size);
      std::uint64_t mapping_data_ptr = 0;
      for (std::uint64_t i = 0; i < 256; ++i)
        for (std::uint64_t j = 0; j < this_superblock_mapping_data[i].size(); ++j)
          superblock_header.m_mapping_data[mapping_data_ptr++] = this_superblock_mapping_data[i][j];
      superblock_header.m_one_bit_count = one_bit_count;
#endif

      // Allocate the variable-size block header.
      superblock_header.m_var_block_headers = std::vector<std::uint8_t>(variable_block_header_size);

      // Fill in variable-size block headers and superblock bitvector.
      std::uint64_t bv_rank = 0;
      sdsl::bit_vector superblock_bv(superblock_bv_size);
      std::vector<std::uint64_t> block_rank(256, 0UL);
      for (std::uint64_t block_id = 0; block_id < blocks_in_this_superblock; ++block_id) {
        block_header_item &block_header = superblock_header.m_block_headers[block_id];
        std::uint64_t block_beg = block_id * block_size;
        std::uint64_t block_end = std::min(block_beg + block_size, this_superblock_size);
        std::uint64_t this_block_size = block_end - block_beg;
        const std::uint8_t *block_ptr = superblock_ptr + block_beg;

        // Fill in the variable-size header of the current block
        // and append the block bitvector to superblock bitvector.
        std::uint64_t ones_count = 0;
        std::uint64_t superblock_bv_offset = block_header.m_bv_offset;
        std::uint64_t variable_block_header_offset = block_header.m_var_size_header_offset;
        encode_block(block_ptr, block_rank, this_block_size, superblock_bv,
            ones_count, superblock_bv_offset,
            superblock_header.m_var_block_headers.data() +
            variable_block_header_offset);

#ifdef ADD_NAVIGATIONAL_BLOCK_HEADER
        block_header.m_bv_rank = bv_rank;
#endif

        // Update block rank.
        bv_rank += ones_count;
        for (std::uint64_t i = 0; i < this_block_size; ++i)
          ++block_rank[block_ptr[i]];
      }

      // Convert the superblock bitvector to final encoding and store.
      superblock_header.m_bitvector = t_bitvector(superblock_bv);
      superblock_header.m_rank_support = t_rank(&superblock_header.m_bitvector);
    }

    void encode_superblock(const std::uint8_t *superblock_ptr, std::uint64_t superblock_id) {
      std::uint64_t superblock_beg = superblock_id * k_superblock_size;
      std::uint64_t superblock_end = std::min(superblock_beg + k_superblock_size, m_size);
      std::uint64_t this_superblock_size = superblock_end - superblock_beg;

      // Store ranks at hyperblock boundary.
      std::uint64_t hyperblock_id = (superblock_id * k_superblock_size) / k_hyperblock_size;
      if (!(superblock_id * k_superblock_size % k_hyperblock_size))
        for (std::uint64_t i = 0; i < 256; ++i)
          m_hyperblock_rank[hyperblock_id * 256 + i] = m_count[i];

      // Store ranks at superblock boundary.
      for (std::uint64_t i = 0; i < 256; ++i)
        m_superblock_rank[superblock_id * 256 + i] =
          m_count[i] - m_hyperblock_rank[hyperblock_id * 256 + i];

      // Update symbol counts.
      for (std::uint64_t i = 0; i < this_superblock_size; ++i)
        ++m_count[(std::uint8_t)superblock_ptr[i]];

      // Compute superblock sigma and mapping from
      // global alphabet to superblock alphabet.
      std::uint64_t this_superblock_sigma = 0;
      for (std::uint64_t i = 0; i < 256; ++i)
        if (m_superblock_rank[superblock_id * 256 + i] + m_hyperblock_rank[hyperblock_id * 256 + i] != m_count[i])
          m_global_mapping[superblock_id * 256 + i] = this_superblock_sigma++;
      m_superblock_headers[superblock_id].m_sigma = this_superblock_sigma - 1;

#ifdef ALLOW_VARIABLE_BLOCK_SIZE
#ifdef FAST_CONSTRUCTION

      // Find the optimal block size. Fast but potentially inoptimal version.
      std::uint64_t best_block_size_log = 0;
      std::uint64_t best_encoding_size = 0;
      std::uint64_t smallest_block_size_log =
        std::max((std::int64_t)0, (std::int64_t)std::min(k_superblock_size_log, (std::uint64_t)16) - 7);
      std::uint64_t smallest_block_size = (1UL << smallest_block_size_log);
      std::uint64_t max_blocks_in_superblock = k_superblock_size / smallest_block_size;

      // Allocate auxiliary arrays used to estimate the
      // size of compressed superblock bitvector.
      std::vector<std::uint64_t> **freq = new std::vector<std::uint64_t>*[max_blocks_in_superblock];
      for (std::uint64_t i = 0; i < max_blocks_in_superblock; ++i)
        freq[i] = new std::vector<std::uint64_t>(256, 0UL);
      std::uint64_t compressed_superblock_bv_size = 0;
      std::uint64_t prev_uncompressed_superblock_bv_size = 0;

      // Try few different block sizes.
      for (std::uint64_t block_size_log = smallest_block_size_log;
          block_size_log <= std::min(k_superblock_size_log, (std::uint64_t)16);
          ++block_size_log) {

        std::uint64_t block_size = (1UL << block_size_log);
        std::uint64_t blocks_in_this_superblock =
          (this_superblock_size + block_size - 1) / block_size;

        // Initialize the encoding size with the space for fixed
        // size block headers and the space for superblock mapping.
        std::uint64_t encoding_size =
          sizeof(block_header_item) * blocks_in_this_superblock +
          this_superblock_sigma * (k_superblock_size / block_size);

        // Compute symbol frequencies for all blocks.
        if (block_size_log == smallest_block_size_log) {

          // Compute the values from scratch.
          for (std::uint64_t block_id = 0; block_id < blocks_in_this_superblock; ++block_id) {
            std::uint64_t block_beg = block_id * block_size;
            std::uint64_t block_end = std::min(block_beg + block_size, this_superblock_size);
            std::uint64_t this_block_size = block_end - block_beg;
            const std::uint8_t *block_ptr = superblock_ptr + block_beg;
            compute_symbol_freq(block_ptr, this_block_size, *freq[block_id]);
          }
        } else {

          // Compute symbol frequency in each block
          // from the frequencies of smaller blocks.
          std::uint64_t prev_blocks_in_this_superblock =
            (this_superblock_size + (block_size / 2) - 1) / (block_size / 2);
          for (std::uint64_t block_id = 0;
              block_id < prev_blocks_in_this_superblock; block_id += 2)
            for (std::uint64_t c = 0; c < 256; ++c)
              (*freq[block_id >> 1])[c] =
                (*freq[block_id])[c] +
                (block_id + 1 < prev_blocks_in_this_superblock ? (*freq[block_id + 1])[c] : 0);
        }

        // Add the space for variable size block headers
        // dependent on block alphabet size.
        for (std::uint64_t block_id = 0; block_id < blocks_in_this_superblock; ++block_id) {
          std::uint64_t block_sigma = 256 - std::count(freq[block_id]->begin(), freq[block_id]->end(), 0UL);
          encoding_size += block_sigma * 4;
#ifdef ADD_NAVIGATIONAL_BLOCK_HEADER
          encoding_size += (block_sigma - 1) * 2;
#endif
        }

        // Compute sizes of block bitvectors for merged blocks.
        std::uint64_t uncompressed_superblock_bv_size = 0;
        for (std::uint64_t block_id = 0; block_id < blocks_in_this_superblock; ++block_id) {

          // Compute Huffman codes.
          std::vector<std::uint64_t> codelen(256);
          compute_huffman_code_lengths(*freq[block_id], codelen);

          // Add the space for variable size block
          // headers dependent on Huffman tree height.
          std::uint64_t max_code_length = *std::max_element(codelen.begin(), codelen.end());
          if (max_code_length > 1)
            encoding_size += (max_code_length - 1) * 3;

          // Compute size of uncompressed block bitvector.
          for (std::uint64_t c = 0; c < 256; ++c)
            uncompressed_superblock_bv_size += (*freq[block_id])[c] * codelen[c];
        }

        if (uncompressed_superblock_bv_size > 0) {

          // We compute the size of headers of compressed superblock
          // bitvector only for the smallest block.
          if (block_size_log == smallest_block_size_log) {
            sdsl::bit_vector block_bv(uncompressed_superblock_bv_size, 0);
            t_bitvector temp_bv(block_bv);
            t_rank temp_rank_support(&temp_bv);
            compressed_superblock_bv_size =
              sdsl::size_in_bytes(temp_bv) +
              sdsl::size_in_bytes(temp_rank_support);
          } else {
            double scaling_factor =
              (double)uncompressed_superblock_bv_size /
              (double)prev_uncompressed_superblock_bv_size;
            compressed_superblock_bv_size =
              (std::uint64_t)((double)compressed_superblock_bv_size * scaling_factor);
          }

          // Add the estimated size of compressed
          // superblock bitvector to total encoding.
          encoding_size += compressed_superblock_bv_size;
        }

        prev_uncompressed_superblock_bv_size = uncompressed_superblock_bv_size;

        // Check if this block size yields better
        // compression than any other tried so far.
        if (block_size_log == smallest_block_size_log ||
            encoding_size < best_encoding_size) {
          best_block_size_log = block_size_log;
          best_encoding_size = encoding_size;
        }
      }

      // Clean up.
      for (std::uint64_t i = 0; i < max_blocks_in_superblock; ++i)
        delete freq[i];
      delete[] freq;

      // Encode blocks inside superblock.
      encode_blocks_in_superblock(superblock_ptr, superblock_id, best_block_size_log);
#else

      // Find the optimal block size. Exact and slow version.
      std::uint64_t best_block_size_log = 0;
      std::uint64_t best_encoding_size = (1UL << 60);
      std::uint64_t smallest_block_size_log =
        std::max((std::int64_t)0, (std::int64_t)std::min(k_superblock_size_log, (std::uint64_t)16) - 7);

      for (std::uint64_t block_size_log = smallest_block_size_log;
          block_size_log <= std::min(k_superblock_size_log, (std::uint64_t)16);
          ++block_size_log) {

        std::uint64_t block_size = (1UL << block_size_log);
        std::uint64_t blocks_in_this_superblock = (this_superblock_size + block_size - 1) / block_size;
        std::uint64_t superblock_beg = superblock_id * k_superblock_size;
        std::uint64_t superblock_end = std::min(superblock_beg + k_superblock_size, m_size);
        std::uint64_t this_superblock_size = superblock_end - superblock_beg;
        std::uint64_t this_superblock_sigma = (std::uint64_t)m_superblock_headers[superblock_id].m_sigma + 1;
        std::uint64_t encoding_size = sizeof(block_header_item) * blocks_in_this_superblock;

#ifdef SPARSE_SUPERBLOCK_MAPPING

        // Create temporary containers for superblock mapping.
        sdsl::bit_vector this_superblock_mapping_bv(
            this_superblock_sigma * (k_superblock_size / block_size), 0);
        std::vector<std::uint8_t> this_superblock_mapping_data[256];
#else
        encoding_size += this_superblock_sigma * (k_superblock_size / block_size);
#endif

        // Compute sizes of variable-size block headers and bitvectors.
        std::uint64_t superblock_bv_size = 0;
        std::uint64_t variable_block_header_size = 0;
        std::vector<std::uint32_t> bv_offset(blocks_in_this_superblock);
        for (std::uint64_t block_id = 0; block_id < blocks_in_this_superblock; ++block_id) {
          std::uint64_t block_beg = block_id * block_size;
          std::uint64_t block_end = std::min(block_beg + block_size, this_superblock_size);
          std::uint64_t this_block_size = block_end - block_beg;
          const std::uint8_t *block_ptr = superblock_ptr + block_beg;

          // Compute global-to-block symbol mapping and basic
          // information about symbol distribution in the block.
          std::uint64_t tree_height = 0, bv_size = 0, sigma = 0;
#ifdef SPARSE_SUPERBLOCK_MAPPING
          std::vector<std::uint16_t> global_to_block_mapping(256, 256);
#endif
          {

            // Compute Huffman code lengths.
            std::vector<std::uint64_t> freq(256), codelen(256);
            compute_symbol_freq(block_ptr, this_block_size, freq);
            compute_huffman_code_lengths(freq, codelen);

            // Sort symbols by frequency.
            std::vector<std::pair<std::uint64_t, std::uint8_t> > sym;
            for (std::uint64_t i = 0; i < 256; ++i)
              if (freq[i] > 0)
                sym.push_back(std::make_pair(codelen[i], (std::uint8_t)i));
            std::sort(sym.begin(), sym.end());

            // Fill in all fields.
            sigma = sym.size();
            tree_height = *std::max_element(codelen.begin(), codelen.end());
            for (std::uint64_t i = 0; i < sym.size(); ++i) {
#ifdef SPARSE_SUPERBLOCK_MAPPING
              global_to_block_mapping[sym[i].second] = (std::uint16_t)i;
#endif
              if (sym.size() > 1)
                bv_size += freq[sym[i].second] * codelen[sym[i].second];
            }
          }

          // Fill in the block header.
          bv_offset[block_id] = (std::uint32_t)superblock_bv_size;

          // Compute sparse representation of superblock mapping.
#ifdef SPARSE_SUPERBLOCK_MAPPING
          for (std::uint64_t i = 0; i < 256; ++i) {
            if (global_to_block_mapping[i] != 256) {
              std::uint8_t superblock_char = m_global_mapping[superblock_id * 256 + i];
              std::uint64_t addr = (std::uint64_t)superblock_char * (k_superblock_size / block_size) + block_id;
              std::uint8_t mapping_val = std::min((std::uint16_t)254, global_to_block_mapping[i]);
              this_superblock_mapping_bv[addr] = 1;
              this_superblock_mapping_data[superblock_char].push_back(mapping_val);
            }
          }
#endif

          // Update the size of superblock bitvector.
          superblock_bv_size += bv_size;

          // Update the total size of variable-size block headers.
          if (tree_height > 1)
            variable_block_header_size += (tree_height - 1) * 3;  // info about each level in the tree
          variable_block_header_size += sigma * 4;                // info about each leaf in the tree
#ifdef ADD_NAVIGATIONAL_BLOCK_HEADER
          variable_block_header_size += (sigma - 1) * 2;          // additional navigational info
#endif
        }

#ifdef SPARSE_SUPERBLOCK_MAPPING
        sdsl::bit_vector_il<> temp_m_mapping_bv(this_superblock_mapping_bv);
        encoding_size += sdsl::size_in_bytes(temp_m_mapping_bv);
        for (std::uint64_t i = 0; i < 256; ++i)
          encoding_size += this_superblock_mapping_data[i].size();
#endif
        encoding_size += variable_block_header_size;

        if (superblock_bv_size > 0) {
          sdsl::bit_vector superblock_bv(superblock_bv_size);
          for (std::uint64_t block_id = 0; block_id < blocks_in_this_superblock; ++block_id) {
            std::uint64_t block_beg = block_id * block_size;
            std::uint64_t block_end = std::min(block_beg + block_size, this_superblock_size);
            std::uint64_t this_block_size = block_end - block_beg;
            const std::uint8_t *block_ptr = superblock_ptr + block_beg;
            compute_block_bv(block_ptr, this_block_size, superblock_bv, (std::uint64_t)bv_offset[block_id]);
          }
          t_bitvector temp_m_bitvector(superblock_bv);
          t_rank temp_m_rank_support(&temp_m_bitvector);
          encoding_size += sdsl::size_in_bytes(temp_m_bitvector);
          encoding_size += sdsl::size_in_bytes(temp_m_rank_support);
        }

        if (encoding_size < best_encoding_size) {
          best_block_size_log = block_size_log;
          best_encoding_size = encoding_size;
        }
      }

      // Encode blocks inside superblock.
      encode_blocks_in_superblock(superblock_ptr, superblock_id, best_block_size_log);
#endif
#else
      encode_blocks_in_superblock(superblock_ptr, superblock_id, k_block_size_log);
#endif
    }

    void init(const std::uint8_t *text, std::uint64_t text_length) {
      m_size = text_length;
      m_count = std::vector<std::uint64_t>(256, 0UL);
      std::uint64_t n_superblocks = (m_size + k_superblock_size - 1) / k_superblock_size;
      std::uint64_t n_hyperblocks = (m_size + k_hyperblock_size - 1) / k_hyperblock_size;

      // Allocate headers.
      m_hyperblock_rank = std::vector<std::uint64_t>(n_hyperblocks * 256);
      m_superblock_rank = std::vector<std::uint32_t>(n_superblocks * 256);
      m_global_mapping = std::vector<std::uint8_t>(n_superblocks * 256, 255);
      m_superblock_headers = std::vector<superblock_header_item>(n_superblocks);

      // Encode superblocks left to right.
      for (std::uint64_t superblock_id = 0; superblock_id < n_superblocks; ++superblock_id) {
        std::uint64_t superblock_beg = superblock_id * k_superblock_size;
        const std::uint8_t *superblock_ptr = text + superblock_beg;
        encode_superblock(superblock_ptr, superblock_id);
      }
    }

    void init(sdsl::int_vector_buffer<(std::uint8_t)8> &text_buf, std::uint64_t text_length) {
      m_size = text_length;
      m_count = std::vector<std::uint64_t>(256, 0UL);
      std::uint64_t n_superblocks = (m_size + k_superblock_size - 1) / k_superblock_size;
      std::uint64_t n_hyperblocks = (m_size + k_hyperblock_size - 1) / k_hyperblock_size;

      // Allocate headers.
      m_hyperblock_rank = std::vector<std::uint64_t>(n_hyperblocks * 256);
      m_superblock_rank = std::vector<std::uint32_t>(n_superblocks * 256);
      m_global_mapping = std::vector<std::uint8_t>(n_superblocks * 256, 255);
      m_superblock_headers = std::vector<superblock_header_item>(n_superblocks);

      // Encode superblocks left to right.
      std::uint8_t *superblock_buf = new std::uint8_t[k_superblock_size];
      for (std::uint64_t superblock_id = 0; superblock_id < n_superblocks; ++superblock_id) {
        std::uint64_t superblock_beg = superblock_id * k_superblock_size;
        std::uint64_t superblock_end = std::min(superblock_beg + k_superblock_size, m_size);
        std::uint64_t this_superblock_size = superblock_end - superblock_beg;
        for (std::uint64_t i = 0; i < this_superblock_size; ++i)
          superblock_buf[i] = text_buf[superblock_beg + i];

        encode_superblock(superblock_buf, superblock_id);
      }
      delete[] superblock_buf;
    }

  public:
    std::uint8_t operator[](std::uint64_t i) const {
      std::uint64_t superblock_id = i / k_superblock_size;
      std::uint64_t superblock_i = i % k_superblock_size;
      const superblock_header_item &superblock_header = m_superblock_headers[superblock_id];

#ifdef ALLOW_VARIABLE_BLOCK_SIZE
      std::uint64_t block_size_log = superblock_header.m_block_size_log;
#else
      std::uint64_t block_size_log = k_block_size_log;
#endif

      std::uint64_t block_size = (1UL << block_size_log);
      std::uint64_t block_i = i & (block_size - 1);
      std::uint64_t this_block_size = std::min(block_size, m_size - (i - block_i));
      std::uint64_t block_id = (superblock_i >> block_size_log);

      const block_header_item &block_header = superblock_header.m_block_headers[block_id];
      std::uint64_t var_size_header_offset = block_header.m_var_size_header_offset;
      std::uint64_t block_tree_height = block_header.m_tree_height;
      const std::uint8_t *variable_block_header_ptr =
        superblock_header.m_var_block_headers.data() + var_size_header_offset;
      const std::uint8_t *copy_variable_block_header_ptr = variable_block_header_ptr;
      const std::uint8_t *variable_block_header_ptr_temp = variable_block_header_ptr;
      if (block_tree_height > 0)
        variable_block_header_ptr_temp += (block_tree_height - 1) * 3;
      const std::uint32_t *variable_block_header_ptr32 =
        (std::uint32_t *)variable_block_header_ptr_temp;

      if (block_tree_height == 0) {
        std::uint8_t c = (variable_block_header_ptr32[0]) & 255;
        return c;
      }

      std::uint64_t code = 0;
      std::uint64_t codelen = 0;

#ifdef ADD_NAVIGATIONAL_BLOCK_HEADER
      std::uint64_t bv_rank = block_header.m_bv_rank;      // rank (in superblock bv) at the beginning of current level
      std::uint64_t bv_offset = block_header.m_bv_offset;  // starting pos (in superblock bv) of bv at current level
      std::uint64_t int_nodes_count = 1;                   // number of internal nodes at current level
      std::uint64_t left_siblings_count = 0;               // number of left siblings (excluding leaves) of current node
      std::uint64_t left_siblings_total_bv_size = 0;       // total size of bv's corresponding to left siblings
                                                           // of current node (excluding leaves)
      std::uint64_t cur_node_bv_size = this_block_size;          // size of bitvector of current node
      std::uint64_t cur_depth_total_bv_size = cur_node_bv_size;  // total size of bitvectors at current depth
      std::uint64_t cur_node_rank = block_i;                     // the rank value refined at each level

      // We maintain second pointer to variable-size block header. It is
      // used to extract total number of 1s in bitvectors corresponding to
      // left siblings of current node (excluding leaves) and also including
      // 1s in the bitvector of current node.
      std::uint64_t block_sigma = (std::uint64_t)block_header.m_sigma + 1;
      variable_block_header_ptr_temp = variable_block_header_ptr;
      if (block_tree_height > 0)
        variable_block_header_ptr_temp += (block_tree_height - 1) * 3;
      variable_block_header_ptr_temp += block_sigma * 4;
      std::uint16_t *second_variable_block_header_ptr =
        (std::uint16_t *)variable_block_header_ptr_temp;

      // Traverse the tree.
      for (std::uint64_t depth = 0; ; ++depth) {

        // Compute the number of 1s in current node.
        std::uint64_t rank1 =
          superblock_header.m_rank_support.rank(bv_offset +
              left_siblings_total_bv_size + cur_node_rank);
        std::uint64_t next_bit =
          superblock_header.m_bitvector[bv_offset + left_siblings_total_bv_size + cur_node_rank];
        std::uint64_t left_siblings_total_ones_count = 0;
        if (left_siblings_count > 0)
          left_siblings_total_ones_count =
            (std::uint64_t)second_variable_block_header_ptr[left_siblings_count - 1];
        rank1 -= bv_rank + left_siblings_total_ones_count;

        // Compute remaining stats about current node.
        std::uint64_t cur_node_one_count =
          second_variable_block_header_ptr[left_siblings_count] -
          left_siblings_total_ones_count;
        std::uint64_t cur_node_zero_count = cur_node_bv_size - cur_node_one_count;
        std::uint64_t rank0 = cur_node_rank - rank1;

        // Update navigational info.
        bv_rank += second_variable_block_header_ptr[int_nodes_count - 1];
        second_variable_block_header_ptr += int_nodes_count;
        left_siblings_count <<= 1;

        // Update rank.
        code <<= 1;
        ++codelen;
        if (next_bit) {
          code |= 1;
          cur_node_rank = rank1;
          cur_node_bv_size = cur_node_one_count;
          ++left_siblings_count;
          left_siblings_total_bv_size += cur_node_zero_count;
        } else {
          cur_node_rank = rank0;
          cur_node_bv_size = cur_node_zero_count;
        }

        // Update navigational info.
        if (depth + 1 < block_tree_height) {

          // Decode leaf count and total size of bitvectors corresponding
          // to internal nodes on the next level of the tree from the
          // variable-size block header.
          std::uint64_t next_level_leaf_count = *variable_block_header_ptr;
          ++variable_block_header_ptr;
          std::uint64_t next_level_total_bv_size =
            (std::uint64_t)(*((std::uint16_t *)variable_block_header_ptr)) + 1;
          variable_block_header_ptr += 2;

          // Update the total size of bitvectors
          // corresponding to left siblings of current node.
          left_siblings_total_bv_size -=
            (cur_depth_total_bv_size - next_level_total_bv_size);

          // Update bitvector offset and total
          // size of bitvectors at current depth.
          bv_offset += cur_depth_total_bv_size;
          cur_depth_total_bv_size = next_level_total_bv_size;

          // Update the number of internal nodes at current level.
          int_nodes_count <<= 1;
          int_nodes_count -= next_level_leaf_count;

          // Update the number of left siblings of current node of exit.
          if (left_siblings_count >= next_level_leaf_count) {
            left_siblings_count -= next_level_leaf_count;
          } else break;
        } else break;
      }
#else
      std::uint64_t cur_depth_bv_offset = block_header.m_bv_offset;
      std::uint64_t cur_depth_bv_total_size = this_block_size;  // total size of bitvectors at current depth
      std::uint64_t cur_node_bv_size = this_block_size;         // size of bitvector of the current node
      std::uint64_t cur_node_rank = block_i;                    // the rank value refined at each level
      std::uint64_t prec_bv_total_size = 0;  // total size of bitvectors to the left of current node at current depth
      bool cached_query = false;
      std::uint64_t cached_query_argument = 0;
      std::uint64_t cached_query_result = 0;

      for (std::uint64_t depth = 0; ; ++depth) {

        // Issue rank query at the beginning, end, and position `block_rank'
        // of the current bitvector. A small optimization is that sometimes
        // the last query issued at a given depth is the same as the first
        // query on the next depth. Thus we cache the last result of rank
        // query.
        std::uint64_t rank_beg = 0;
        if (cached_query && cached_query_argument == cur_depth_bv_offset + prec_bv_total_size)
          rank_beg = cached_query_result;
        else rank_beg = superblock_header.m_rank_support.rank(cur_depth_bv_offset + prec_bv_total_size);
        std::uint64_t rank_mid = superblock_header.m_rank_support.rank(
            cur_depth_bv_offset + prec_bv_total_size + cur_node_rank);
        std::uint64_t rank_end = superblock_header.m_rank_support.rank(
            cur_depth_bv_offset + prec_bv_total_size + cur_node_bv_size);
        std::uint64_t next_bit =
          superblock_header.m_bitvector[cur_depth_bv_offset + prec_bv_total_size + cur_node_rank];
        cached_query_argument = cur_depth_bv_offset + prec_bv_total_size + cur_node_bv_size;
        cached_query_result = rank_end;
        cached_query = true;

        // The number of 0/1 bits in the current internal node.
        std::uint64_t cur_node_bv_one_count = rank_end - rank_beg;
        std::uint64_t cur_node_bv_zero_count = cur_node_bv_size - cur_node_bv_one_count;

        // The number of 0/1 bits up to queried position.
        std::uint64_t rank1 = rank_mid - rank_beg;
        std::uint64_t rank0 = cur_node_rank - rank1;

        // Update rank.
        ++codelen;
        code <<= 1;
        if (next_bit) {
          code |= 1;
          cur_node_rank = rank1;
          cur_node_bv_size = cur_node_bv_one_count;
          prec_bv_total_size += cur_node_bv_zero_count;
        } else {
          cur_node_rank = rank0;
          cur_node_bv_size = cur_node_bv_zero_count;
        }

        // Update navigational info.
        if (depth + 1 < block_tree_height) {
          ++variable_block_header_ptr;
          std::uint64_t next_level_bv_total_size =
            (std::uint64_t)(*((std::uint16_t *)variable_block_header_ptr)) + 1;
          variable_block_header_ptr += 2;

          if (prec_bv_total_size >= cur_depth_bv_total_size - next_level_bv_total_size)
            prec_bv_total_size -= (cur_depth_bv_total_size - next_level_bv_total_size);
          else break;

          cur_depth_bv_offset += cur_depth_bv_total_size;
          cur_depth_bv_total_size = next_level_bv_total_size;
        } else break;
      }
#endif

      std::uint64_t block_c = compute_symbol_from_block_header(
          copy_variable_block_header_ptr, code, codelen);
      std::uint8_t c = (variable_block_header_ptr32[block_c] & 255);

      return c;
    }

    std::pair<std::uint64_t, std::uint8_t>
    inverse_select(std::uint64_t i) const {
      std::uint64_t hyperblock_id = i / k_hyperblock_size;
      std::uint64_t superblock_id = i / k_superblock_size;
      std::uint64_t superblock_i = i % k_superblock_size;
      const superblock_header_item &superblock_header = m_superblock_headers[superblock_id];

#ifdef ALLOW_VARIABLE_BLOCK_SIZE
      std::uint64_t block_size_log = superblock_header.m_block_size_log;
#else
      std::uint64_t block_size_log = k_block_size_log;
#endif

      std::uint64_t block_size = (1UL << block_size_log);
      std::uint64_t block_i = i & (block_size - 1);
      std::uint64_t this_block_size = std::min(block_size, m_size - (i - block_i));
      std::uint64_t block_id = (superblock_i >> block_size_log);

      const block_header_item &block_header = superblock_header.m_block_headers[block_id];
      std::uint64_t var_size_header_offset = block_header.m_var_size_header_offset;
      std::uint64_t block_tree_height = block_header.m_tree_height;
      const std::uint8_t *variable_block_header_ptr =
        superblock_header.m_var_block_headers.data() + var_size_header_offset;
      const std::uint8_t *copy_variable_block_header_ptr = variable_block_header_ptr;
      const std::uint8_t *variable_block_header_ptr_temp = variable_block_header_ptr;
      if (block_tree_height > 0)
        variable_block_header_ptr_temp += (block_tree_height - 1) * 3;
      const std::uint32_t *variable_block_header_ptr32 =
        (std::uint32_t *)variable_block_header_ptr_temp;

      if (block_tree_height == 0) {
        std::uint64_t block_c = 0;
        std::uint8_t c = (variable_block_header_ptr32[block_c]) & 255;

        if (i == 0)
          return std::make_pair((std::uint64_t)0, c);
        else {
          std::uint64_t cur_node_rank = block_i;
          std::uint64_t rank_at_block_boundary = (variable_block_header_ptr32[block_c] >> 8);
          std::uint64_t rank_at_superblock_boundary = m_superblock_rank[superblock_id * 256 + c];
          std::uint64_t rank_at_hyperblock_boundary = m_hyperblock_rank[hyperblock_id * 256 + c];

          std::uint64_t ret_rank =
            rank_at_hyperblock_boundary +
            rank_at_superblock_boundary +
            rank_at_block_boundary +
            cur_node_rank;

          return std::make_pair(ret_rank, c);
        }
      }

      std::uint64_t code = 0;
      std::uint64_t codelen = 0;

#ifdef ADD_NAVIGATIONAL_BLOCK_HEADER
      std::uint64_t bv_rank = block_header.m_bv_rank;      // rank (in superblock bv) at the beginning of current level
      std::uint64_t bv_offset = block_header.m_bv_offset;  // starting pos (in superblock bv) of bv at current level
      std::uint64_t int_nodes_count = 1;                   // number of internal nodes at current level
      std::uint64_t left_siblings_count = 0;               // number of left siblings (excluding leaves) of current node
      std::uint64_t left_siblings_total_bv_size = 0;       // total size of bv's corresponding to left siblings
                                                           // of current node (excluding leaves)
      std::uint64_t cur_node_bv_size = this_block_size;          // size of bitvector of current node
      std::uint64_t cur_depth_total_bv_size = cur_node_bv_size;  // total size of bitvectors at current depth
      std::uint64_t cur_node_rank = block_i;                     // the rank value refined at each level

      // We maintain second pointer to variable-size block header. It is
      // used to extract total number of 1s in bitvectors corresponding to
      // left siblings of current node (excluding leaves) and also including
      // 1s in the bitvector of current node.
      std::uint64_t block_sigma = (std::uint64_t)block_header.m_sigma + 1;
      variable_block_header_ptr_temp = variable_block_header_ptr;
      if (block_tree_height > 0)
        variable_block_header_ptr_temp += (block_tree_height - 1) * 3;
      variable_block_header_ptr_temp += block_sigma * 4;
      std::uint16_t *second_variable_block_header_ptr =
        (std::uint16_t *)variable_block_header_ptr_temp;

      // Traverse the tree.
      for (std::uint64_t depth = 0; ; ++depth) {

        // Compute the number of 1s in current node.
        std::uint64_t rank1 =
          superblock_header.m_rank_support.rank(bv_offset +
              left_siblings_total_bv_size + cur_node_rank);
        std::uint64_t next_bit =
          superblock_header.m_bitvector[bv_offset + left_siblings_total_bv_size + cur_node_rank];
        std::uint64_t left_siblings_total_ones_count = 0;
        if (left_siblings_count > 0)
          left_siblings_total_ones_count =
            (std::uint64_t)second_variable_block_header_ptr[left_siblings_count - 1];
        rank1 -= bv_rank + left_siblings_total_ones_count;

        // Compute remaining stats about current node.
        std::uint64_t cur_node_one_count =
          second_variable_block_header_ptr[left_siblings_count] -
          left_siblings_total_ones_count;
        std::uint64_t cur_node_zero_count = cur_node_bv_size - cur_node_one_count;
        std::uint64_t rank0 = cur_node_rank - rank1;

        // Update navigational info.
        bv_rank += second_variable_block_header_ptr[int_nodes_count - 1];
        second_variable_block_header_ptr += int_nodes_count;
        left_siblings_count <<= 1;

        // Update rank.
        code <<= 1;
        ++codelen;
        if (next_bit) {
          code |= 1;
          cur_node_rank = rank1;
          cur_node_bv_size = cur_node_one_count;
          ++left_siblings_count;
          left_siblings_total_bv_size += cur_node_zero_count;
        } else {
          cur_node_rank = rank0;
          cur_node_bv_size = cur_node_zero_count;
        }

        // Update navigational info.
        if (depth + 1 < block_tree_height) {

          // Decode leaf count and total size of bitvectors corresponding
          // to internal nodes on the next level of the tree from the
          // variable-size block header.
          std::uint64_t next_level_leaf_count = *variable_block_header_ptr;
          ++variable_block_header_ptr;
          std::uint64_t next_level_total_bv_size =
            (std::uint64_t)(*((std::uint16_t *)variable_block_header_ptr)) + 1;
          variable_block_header_ptr += 2;

          // Update the total size of bitvectors
          // corresponding to left siblings of current node.
          left_siblings_total_bv_size -=
            (cur_depth_total_bv_size - next_level_total_bv_size);

          // Update bitvector offset and total
          // size of bitvectors at current depth.
          bv_offset += cur_depth_total_bv_size;
          cur_depth_total_bv_size = next_level_total_bv_size;

          // Update the number of internal nodes at current level.
          int_nodes_count <<= 1;
          int_nodes_count -= next_level_leaf_count;

          // Update the number of left siblings of current node of exit.
          if (left_siblings_count >= next_level_leaf_count) {
            left_siblings_count -= next_level_leaf_count;
          } else break;
        } else break;
      }
#else
      std::uint64_t cur_depth_bv_offset = block_header.m_bv_offset;
      std::uint64_t cur_depth_bv_total_size = this_block_size;  // total size of bitvectors at current depth
      std::uint64_t cur_node_bv_size = this_block_size;         // size of bitvector of the current node
      std::uint64_t cur_node_rank = block_i;                    // the rank value refined at each level
      std::uint64_t prec_bv_total_size = 0;  // total size of bitvectors to the left of current node at current depth
      bool cached_query = false;
      std::uint64_t cached_query_argument = 0;
      std::uint64_t cached_query_result = 0;

      for (std::uint64_t depth = 0; ; ++depth) {

        // Issue rank query at the beginning, end, and position `block_rank'
        // of the current bitvector. A small optimization is that sometimes
        // the last query issued at a given depth is the same as the first
        // query on the next depth. Thus we cache the last result of rank
        // query.
        std::uint64_t rank_beg = 0;
        if (cached_query && cached_query_argument == cur_depth_bv_offset + prec_bv_total_size)
          rank_beg = cached_query_result;
        else rank_beg = superblock_header.m_rank_support.rank(cur_depth_bv_offset + prec_bv_total_size);
        std::uint64_t rank_mid = superblock_header.m_rank_support.rank(
            cur_depth_bv_offset + prec_bv_total_size + cur_node_rank);
        std::uint64_t rank_end = superblock_header.m_rank_support.rank(
            cur_depth_bv_offset + prec_bv_total_size + cur_node_bv_size);
        std::uint64_t next_bit =
          superblock_header.m_bitvector[cur_depth_bv_offset + prec_bv_total_size + cur_node_rank];
        cached_query_argument = cur_depth_bv_offset + prec_bv_total_size + cur_node_bv_size;
        cached_query_result = rank_end;
        cached_query = true;

        // The number of 0/1 bits in the current internal node.
        std::uint64_t cur_node_bv_one_count = rank_end - rank_beg;
        std::uint64_t cur_node_bv_zero_count = cur_node_bv_size - cur_node_bv_one_count;

        // The number of 0/1 bits up to queried position.
        std::uint64_t rank1 = rank_mid - rank_beg;
        std::uint64_t rank0 = cur_node_rank - rank1;

        // Update rank.
        ++codelen;
        code <<= 1;
        if (next_bit) {
          code |= 1;
          cur_node_rank = rank1;
          cur_node_bv_size = cur_node_bv_one_count;
          prec_bv_total_size += cur_node_bv_zero_count;
        } else {
          cur_node_rank = rank0;
          cur_node_bv_size = cur_node_bv_zero_count;
        }

        // Update navigational info.
        if (depth + 1 < block_tree_height) {
          ++variable_block_header_ptr;
          std::uint64_t next_level_bv_total_size =
            (std::uint64_t)(*((std::uint16_t *)variable_block_header_ptr)) + 1;
          variable_block_header_ptr += 2;

          if (prec_bv_total_size >= cur_depth_bv_total_size - next_level_bv_total_size)
            prec_bv_total_size -= (cur_depth_bv_total_size - next_level_bv_total_size);
          else break;

          cur_depth_bv_offset += cur_depth_bv_total_size;
          cur_depth_bv_total_size = next_level_bv_total_size;
        } else break;
      }
#endif

      std::uint64_t block_c = compute_symbol_from_block_header(
          copy_variable_block_header_ptr, code, codelen);
      std::uint8_t c = (variable_block_header_ptr32[block_c] & 255);

      if (i == 0)
        return std::make_pair((std::uint64_t)0, c);
      else {
        std::uint64_t rank_at_block_boundary = (variable_block_header_ptr32[block_c] >> 8);
        std::uint64_t rank_at_superblock_boundary = m_superblock_rank[superblock_id * 256 + c];
        std::uint64_t rank_at_hyperblock_boundary = m_hyperblock_rank[hyperblock_id * 256 + c];

        std::uint64_t ret_rank =
          rank_at_hyperblock_boundary +
          rank_at_superblock_boundary +
          rank_at_block_boundary +
          cur_node_rank;

        return std::make_pair(ret_rank, c);
      }
    }

    // We assume 0 <= i <= m_size.
    std::uint64_t rank(std::uint64_t i, std::uint8_t c) const {
      if (i == 0) return 0;
      else if (i == m_size) return m_count[c];

      std::uint64_t hyperblock_id = i / k_hyperblock_size;
      std::uint64_t superblock_id = i / k_superblock_size;
      std::uint64_t superblock_i = i % k_superblock_size;
      std::uint8_t superblock_c = m_global_mapping[superblock_id * 256 + c];
      const superblock_header_item &superblock_header = m_superblock_headers[superblock_id];
      std::uint64_t superblock_sigma = (std::uint64_t)superblock_header.m_sigma + 1;

#ifdef ALLOW_VARIABLE_BLOCK_SIZE
      std::uint64_t block_size_log = superblock_header.m_block_size_log;
#else
      std::uint64_t block_size_log = k_block_size_log;
#endif

      std::uint64_t block_size = (1UL << block_size_log);
      std::uint64_t blocks_in_superblock_log = (k_superblock_size_log - block_size_log);
      std::uint64_t block_i = i & (block_size - 1);
      std::uint64_t this_block_size = std::min(block_size, m_size - (i - block_i));
      std::uint64_t block_id = (superblock_i >> block_size_log);
      std::uint64_t rank_at_superblock_boundary = m_superblock_rank[superblock_id * 256 + c];
      std::uint64_t rank_at_hyperblock_boundary = m_hyperblock_rank[hyperblock_id * 256 + c];

      if (superblock_c >= superblock_sigma)  // special case: c does not occur in the superblock
        return
          rank_at_hyperblock_boundary +
          rank_at_superblock_boundary;

#ifdef SPARSE_SUPERBLOCK_MAPPING
      std::uint64_t addr = ((std::uint64_t)superblock_c << blocks_in_superblock_log) + block_id;
      bool mapping_bit = superblock_header.m_mapping_bv[addr];
      std::uint64_t mapping_rank = superblock_header.m_mapping_bv_rank_support.rank(addr);
      std::uint8_t block_c = mapping_bit ? superblock_header.m_mapping_data[mapping_rank] : 255;
#else
      const std::uint8_t *superblock_mapping_ptr = superblock_header.m_mapping.data();
      std::uint8_t block_c =
        superblock_mapping_ptr[((std::uint64_t)superblock_c << blocks_in_superblock_log) + block_id];
#endif

      if (block_c == 255) {  // special case: c does not occur in the block
#ifdef SPARSE_SUPERBLOCK_MAPPING
        std::uint64_t select_ret = 0;
        if (mapping_rank != superblock_header.m_one_bit_count)
          select_ret = superblock_header.m_mapping_bv_select_support.select(mapping_rank + 1);
        if (mapping_rank == superblock_header.m_one_bit_count ||
            select_ret >= ((superblock_c + 1) << blocks_in_superblock_log)) {
#else

        // Find the closest block to the right in which c occurs.
        ++block_id;
        std::uint64_t blocks_in_superblock = (1UL << blocks_in_superblock_log);
        while (block_id < blocks_in_superblock &&
            superblock_mapping_ptr[((std::uint64_t)superblock_c << blocks_in_superblock_log) + block_id] == 255)
          ++block_id;
        if (block_id == blocks_in_superblock) {
#endif

          // Return the answer from superblock header or global counts.
          if ((superblock_id + 1) * k_superblock_size >= m_size) return m_count[c];
          else return rank_at_hyperblock_boundary + m_superblock_rank[(superblock_id + 1) * 256 + c];
        } else {
#ifdef SPARSE_SUPERBLOCK_MAPPING
          block_c = superblock_header.m_mapping_data[mapping_rank];
          block_id = superblock_header.m_mapping_bv_select_support.select(mapping_rank + 1) -
            ((std::uint64_t)superblock_c << blocks_in_superblock_log);
#else
          block_c =
            superblock_mapping_ptr[((std::uint64_t)superblock_c << blocks_in_superblock_log) + block_id];
#endif

          // Return the rank value from block header.
          const block_header_item &block_header = superblock_header.m_block_headers[block_id];
          std::uint64_t var_size_header_offset = block_header.m_var_size_header_offset;
          std::uint64_t block_tree_height = block_header.m_tree_height;
          const std::uint8_t *variable_block_header_ptr =
            superblock_header.m_var_block_headers.data() + var_size_header_offset;
          if (block_tree_height > 0)
            variable_block_header_ptr += (block_tree_height - 1) * 3;
          std::uint32_t *variable_block_header_ptr32 = (std::uint32_t *)variable_block_header_ptr;
          if ((variable_block_header_ptr32[block_c] & 0xff) != c)  // special case: block_c == 254
            ++block_c;
          std::uint64_t rank_at_block_boundary = (variable_block_header_ptr32[block_c] >> 8);
          return
            rank_at_hyperblock_boundary +
            rank_at_superblock_boundary +
            rank_at_block_boundary;
        }
      }

      // Compute rank at block boundary.
      const block_header_item &block_header = superblock_header.m_block_headers[block_id];
      std::uint64_t var_size_header_offset = block_header.m_var_size_header_offset;
      std::uint64_t block_tree_height = block_header.m_tree_height;
      const std::uint8_t *variable_block_header_ptr =
        superblock_header.m_var_block_headers.data() + var_size_header_offset;
      const std::uint8_t *variable_block_header_ptr_temp = variable_block_header_ptr;
      if (block_tree_height > 0)
        variable_block_header_ptr_temp += (block_tree_height - 1) * 3;
      const std::uint32_t *variable_block_header_ptr32 = (std::uint32_t *)variable_block_header_ptr_temp;
      if ((variable_block_header_ptr32[block_c] & 255) != c)  // special case: block_c == 254
        ++block_c;
      std::uint64_t rank_at_block_boundary = (variable_block_header_ptr32[block_c] >> 8);

      // Answer rank query inside block.
      if (block_tree_height == 0) // special case: block was a run of single symbol
        return
          rank_at_hyperblock_boundary +
          rank_at_superblock_boundary +
          rank_at_block_boundary +
          block_i;

      std::uint64_t code = 0;
      std::uint64_t codelen = 0;
      restore_code_from_block_header(block_c,
          variable_block_header_ptr, block_tree_height, code, codelen);

#ifdef ADD_NAVIGATIONAL_BLOCK_HEADER
      std::uint64_t bv_rank = block_header.m_bv_rank;        // rank (in superblock bv) at the beginning of current level
      std::uint64_t bv_offset = block_header.m_bv_offset;    // starting pos (in superblock bv) of bv at current level
      std::uint64_t int_nodes_count = 1;                     // number of internal nodes at current level
      std::uint64_t left_siblings_count = 0;                 // number of left siblings (excluding leaves) of current node
      std::uint64_t left_siblings_total_bv_size = 0;         // total size of bv's corresponding to left siblings of
                                                             // current node (excluding leaves)
      std::uint64_t cur_node_bv_size = this_block_size;          // size of bitvector of current node
      std::uint64_t cur_depth_total_bv_size = cur_node_bv_size;  // total size of bitvectors at current depth
      std::uint64_t cur_node_rank = block_i;                     // the rank value refined at each level

      // We maintain second pointer to variable-size block header. It is
      // used to extract total number of 1s in bitvectors corresponding to
      // left siblings of current node (excluding leaves) and also including
      // 1s in the bitvector of current node.
      std::uint64_t block_sigma = (std::uint64_t)block_header.m_sigma + 1;
      variable_block_header_ptr_temp = variable_block_header_ptr;
      if (block_tree_height > 0)
        variable_block_header_ptr_temp += (block_tree_height - 1) * 3;
      variable_block_header_ptr_temp += block_sigma * 4;
      std::uint16_t *second_variable_block_header_ptr =
        (std::uint16_t *)variable_block_header_ptr_temp;

      // Traverse the tree.
      for (std::uint64_t depth = 0; depth < codelen; ++depth) {

        // Compute the number of 1s in current node.
        std::uint64_t rank1 =
          superblock_header.m_rank_support.rank(bv_offset +
              left_siblings_total_bv_size + cur_node_rank);
        std::uint64_t left_siblings_total_ones_count = 0;
        if (left_siblings_count > 0)
          left_siblings_total_ones_count =
            (std::uint64_t)second_variable_block_header_ptr[left_siblings_count - 1];
        rank1 -= bv_rank + left_siblings_total_ones_count;

        // Compute remaining stats about current node.
        std::uint64_t cur_node_one_count =
          second_variable_block_header_ptr[left_siblings_count] -
          left_siblings_total_ones_count;
        std::uint64_t cur_node_zero_count = cur_node_bv_size - cur_node_one_count;
        std::uint64_t rank0 = cur_node_rank - rank1;

        // Update navigational info.
        bv_rank += second_variable_block_header_ptr[int_nodes_count - 1];
        second_variable_block_header_ptr += int_nodes_count;
        left_siblings_count <<= 1;

        // Update rank.
        std::uint64_t next_bit = (code & (1UL << (codelen - depth - 1)));
        if (next_bit) {
          cur_node_rank = rank1;
          cur_node_bv_size = cur_node_one_count;
          ++left_siblings_count;
          left_siblings_total_bv_size += cur_node_zero_count;
        } else {
          cur_node_rank = rank0;
          cur_node_bv_size = cur_node_zero_count;
        }

        // Update navigational info.
        if (depth + 1 != codelen) {

          // Decode leaf count and total size of bitvectors corresponding
          // to internal nodes on the next level of the tree from the
          // variable-size block header.
          std::uint64_t next_level_leaf_count = *variable_block_header_ptr;
          ++variable_block_header_ptr;
          std::uint64_t next_level_total_bv_size =
            (std::uint64_t)(*((std::uint16_t *)variable_block_header_ptr)) + 1;
          variable_block_header_ptr += 2;

          // Update the total size of bitvectors
          // corresponding to left siblings of current node.
          left_siblings_total_bv_size -=
            (cur_depth_total_bv_size - next_level_total_bv_size);

          // Update bitvector offset and total
          // size of bitvectors at current depth.
          bv_offset += cur_depth_total_bv_size;
          cur_depth_total_bv_size = next_level_total_bv_size;

          // Update the number of internal nodes at current level.
          int_nodes_count <<= 1;
          int_nodes_count -= next_level_leaf_count;

          // Update the number of left siblings of current node.
          left_siblings_count -= next_level_leaf_count;
        }
      }
#else
      std::uint64_t cur_depth_bv_offset = block_header.m_bv_offset;
      std::uint64_t cur_depth_bv_total_size = this_block_size;  // total size of bitvectors at current depth
      std::uint64_t cur_node_bv_size = this_block_size;         // size of bitvector of the current node
      std::uint64_t cur_node_rank = block_i;                    // the rank value refined at each level
      std::uint64_t prec_bv_total_size = 0;                     // total size of bitvectors to the left
                                                                // of current node at current depth
      bool cached_query = false;
      std::uint64_t cached_query_argument = 0;
      std::uint64_t cached_query_result = 0;
      for (std::uint64_t depth = 0; depth < codelen; ++depth) {

        // Issue rank query at the beginning, end, and position `block_rank'
        // of the current bitvector. A small optimization is that sometimes
        // the last query issued at a given depth is the same as the first
        // query on the next depth. Thus we cache the last result of rank
        // query.
        std::uint64_t rank_beg = 0;
        if (cached_query && cached_query_argument == cur_depth_bv_offset + prec_bv_total_size)
          rank_beg = cached_query_result;
        else rank_beg = superblock_header.m_rank_support.rank(
            cur_depth_bv_offset + prec_bv_total_size);
        std::uint64_t rank_mid = superblock_header.m_rank_support.rank(
            cur_depth_bv_offset + prec_bv_total_size + cur_node_rank);
        std::uint64_t rank_end = superblock_header.m_rank_support.rank(
            cur_depth_bv_offset + prec_bv_total_size + cur_node_bv_size);
        cached_query_argument = cur_depth_bv_offset +
          prec_bv_total_size + cur_node_bv_size;
        cached_query_result = rank_end;
        cached_query = true;

        // The number of 0/1 bits in the current internal node.
        std::uint64_t cur_node_bv_one_count = rank_end - rank_beg;
        std::uint64_t cur_node_bv_zero_count =
          cur_node_bv_size - cur_node_bv_one_count;

        // The number of 0/1 bits up to queried position.
        std::uint64_t rank1 = rank_mid - rank_beg;
        std::uint64_t rank0 = cur_node_rank - rank1;

        // Update rank.
        if (code & (1UL << (codelen - depth - 1))) {
          cur_node_rank = rank1;
          cur_node_bv_size = cur_node_bv_one_count;
          prec_bv_total_size += cur_node_bv_zero_count;
        } else {
          cur_node_rank = rank0;
          cur_node_bv_size = cur_node_bv_zero_count;
        }

        // Update navigational info.
        if (depth + 1 != codelen) {
          ++variable_block_header_ptr;
          std::uint64_t next_level_bv_total_size =
            (std::uint64_t)(*((std::uint16_t *)variable_block_header_ptr)) + 1;
          variable_block_header_ptr += 2;
          prec_bv_total_size -=
            (cur_depth_bv_total_size - next_level_bv_total_size);
          cur_depth_bv_offset += cur_depth_bv_total_size;
          cur_depth_bv_total_size = next_level_bv_total_size;
        }
      }
#endif

      return
        rank_at_hyperblock_boundary +
        rank_at_superblock_boundary +
        rank_at_block_boundary +
        cur_node_rank;
    }

    // Swap method
    void swap(wt_fbb& tree) {
      if (this != &tree) {
        std::swap(m_size, tree.m_size);
        std::swap(m_hyperblock_rank, tree.m_hyperblock_rank);
        std::swap(m_superblock_rank, tree.m_superblock_rank);
        std::swap(m_global_mapping, tree.m_global_mapping);
        std::swap(m_superblock_headers, tree.m_superblock_headers);
        std::swap(m_count, tree.m_count);
      }
    }

    // Assignment operator
    wt_fbb& operator=(const wt_fbb& tree) {
      if (this != &tree)
        copy(tree);
      return *this;
    }

    // Move assignment operator
    wt_fbb& operator=(wt_fbb&& tree) {
      swap(tree);
      return *this;
    }

    // Returns the size of the original sequence
    std::uint64_t size() const {
      return m_size;
    }

    // Serializes the data structure into a stream
    std::uint64_t serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const {
      sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
      std::uint64_t written_bytes = 0;
      written_bytes += sdsl::serialize(m_size, out, child, "size");
      written_bytes += sdsl::serialize(m_count, out, child, "count");
      written_bytes += sdsl::serialize(m_hyperblock_rank, out, child, "hyperblock_rank");
      written_bytes += sdsl::serialize(m_superblock_rank, out, child, "superblock_rank");
      written_bytes += sdsl::serialize(m_global_mapping, out, child, "global_mapping");
      written_bytes += sdsl::serialize(m_superblock_headers, out, child, "superblock_headers");
      sdsl::structure_tree::add_size(child, written_bytes);
      return written_bytes;
    }

    // Load the data structure from a stream and set the supported vector
    void load(std::istream& in) {
      sdsl::load(m_size, in);
      sdsl::load(m_count, in);
      sdsl::load(m_hyperblock_rank, in);
      sdsl::load(m_superblock_rank, in);
      sdsl::load(m_global_mapping, in);
      sdsl::load(m_superblock_headers, in);
    }
};

#ifdef ALLOW_VARIABLE_BLOCK_SIZE
template<class t_bitvector, class t_rank, std::uint64_t t_sbs_log>
  const std::uint64_t wt_fbb<t_bitvector, t_rank, t_sbs_log>::k_superblock_size_log = t_sbs_log;
template<class t_bitvector, class t_rank,std::uint64_t t_sbs_log>
  const std::uint64_t wt_fbb<t_bitvector, t_rank, t_sbs_log>::k_superblock_size = (1UL << t_sbs_log);
template<class t_bitvector, class t_rank, std::uint64_t t_sbs_log>
  const std::uint64_t wt_fbb<t_bitvector, t_rank, t_sbs_log>::k_hyperblock_size = (1UL << 32);
#else
template<class t_bitvector, class t_rank, std::uint64_t t_bs_log, std::uint64_t t_sbs_log>
  const std::uint64_t wt_fbb<t_bitvector, t_rank, t_bs_log, t_sbs_log>::k_block_size_log = t_bs_log;
template<class t_bitvector, class t_rank, std::uint64_t t_bs_log, std::uint64_t t_sbs_log>
  const std::uint64_t wt_fbb<t_bitvector, t_rank, t_bs_log, t_sbs_log>::k_block_size = (1UL << t_bs_log);
template<class t_bitvector, class t_rank, std::uint64_t t_bs_log, std::uint64_t t_sbs_log>
  const std::uint64_t wt_fbb<t_bitvector, t_rank, t_bs_log, t_sbs_log>::k_superblock_size_log = t_sbs_log;
template<class t_bitvector, class t_rank, std::uint64_t t_bs_log, std::uint64_t t_sbs_log>
  const std::uint64_t wt_fbb<t_bitvector, t_rank, t_bs_log, t_sbs_log>::k_superblock_size = (1UL << t_sbs_log);
template<class t_bitvector, class t_rank, std::uint64_t t_bs_log, std::uint64_t t_sbs_log>
  const std::uint64_t wt_fbb<t_bitvector, t_rank, t_bs_log, t_sbs_log>::k_hyperblock_size = (1UL << 32);
#endif

#endif  // __WT_FBB_HPP_INCLUDED
