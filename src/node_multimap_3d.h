///////////////////////////////////////////////////////////////////////////////
///
///	\file    node_multimap_3d.h
///	\author  Paul Ullrich
///	\version March 25, 2019
///
///	<remarks>
///		Copyright 2019 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#ifndef _NODEMULTIMAP3D_H_
#define _NODEMULTIMAP3D_H_

#include <iostream>
#include <array>
#include <unordered_set>

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class for hashing Nodes with floating point tolerance.
///	</summary>
template <
	class NodeType,
	class DataType
> class node_multimap_3d {
public:
	typedef std::pair<NodeType, DataType> value_type;

	//typedef std::array<int,3> bin_index_type;

	typedef size_t bin_index_type;

	typedef std::vector<value_type> bin_contents_type;

	typedef std::map<bin_index_type, bin_contents_type> node_map_type;

	typedef typename node_map_type::iterator node_map_type_iterator;

	static const size_t ix_not_found = (-1);

public:
	///	<summary>
	///		A const_iterator to the contents of this node_map_3d.
	///	</summary>
	class const_iterator {
		public:
			const node_map_type * parent_map;
			typename node_map_type::const_iterator iter_node_map;
			typename bin_contents_type::const_iterator iter_bin_contents;

		public:
			const_iterator(
				const node_map_type * _parent_map
			) :
				parent_map(_parent_map)
			{
			}

			const_iterator & increment() {
				//if (parent_map == NULL) {
				//	_EXCEPTION();
				//}
				iter_bin_contents++;
				if (iter_bin_contents == iter_node_map->second.end()) {
					iter_node_map++;
					if (iter_node_map != parent_map->end()) {
						iter_bin_contents = iter_node_map->second.begin();
					}
				}
				return *this;
			}

			const_iterator & operator++() {
				return increment();
			}

			const_iterator & operator++(int) {
				return increment();
			}

			bool operator==(const const_iterator & iter) {
				if (parent_map == iter.parent_map) {
					if (iter_node_map == iter.iter_node_map) {
						if (iter_node_map == parent_map->end()) {
							return true;
						} else {
							if (iter_bin_contents == iter.iter_bin_contents) {
								return true;
							}
						}
					}
				}
				return false;
			}

			bool operator!=(const const_iterator & iter) {
				return !((*this)==iter);
			}

			const value_type & operator*() {
				return *iter_bin_contents;
			}

			const value_type * operator->() {
				return &(*iter_bin_contents);
			}
	};

private:
	///	<summary>
	///		Check if a given bin holds the given Node.
	///	</summary>
	size_t _bin_contains(
		const bin_contents_type & vec,
		const NodeType & node
	) const {
		for (size_t i = 0; i < vec.size(); i++) {
			if (vec[i].first == node) {
				return i;
			}
		}
		return (ix_not_found);
	}

	///	<summary>
	///		Check if a given bin range holds the given Node.
	///	</summary>
	const_iterator _find_in_bin_range(
		const std::array<int,3> & arrBegin,
		const std::array<int,3> & arrEnd,
		const NodeType & node
	) const {
		std::array<int,3> arrCurrent;
		for (int i = arrBegin[0]; i < arrEnd[0]; i++) {
		for (int j = arrBegin[1]; j < arrEnd[1]; j++) {
		for (int k = arrBegin[2]; k < arrEnd[2]; k++) {
			arrCurrent[0] = i;
			arrCurrent[1] = j;
			arrCurrent[2] = k;

			typename node_map_type::const_iterator iter = m_map.find(hasher(arrCurrent));
			if (iter != m_map.end()) {
				size_t ix = _bin_contains(iter->second, node);
				if (ix != ix_not_found) {
					const_iterator iter_out(&m_map);
					iter_out.iter_node_map = iter;
					iter_out.iter_bin_contents = iter->second.begin() + ix;
					return iter_out;
				}
			}
		}
		}
		}

		return end();
	}

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	node_multimap_3d(
		double tolerance = 1.0e-12,
		double bin_width = 1.0e-1
	) :
		m_tolerance(tolerance),
		m_bin_width(bin_width)
	{ }

	///	<summary>
	///		Hash from array to index.
	///	</summary>
	bin_index_type hasher(
		const std::array<int,3> & arr
	) const {
		return std::hash<int>{}(arr[0] * 18397 + arr[1] * 20483 + arr[2] * 29303);
	}

	///	<summary>
	///		Begin const_iterator of this class.
	///	</summary>
	const_iterator begin() const {
		const_iterator iter(&m_map);
		iter.iter_node_map = m_map.begin();
		if (m_map.begin() != m_map.end()) {
			iter.iter_bin_contents = m_map.begin()->second.begin();
		}
		return iter;
	}

	///	<summary>
	///		End const_iterator of this class.
	///	</summary>
	const_iterator end() const {
		const_iterator iter(&m_map);
		iter.iter_node_map = m_map.end();
		return iter;
	}

	///	<summary>
	///		Number of bins in this object.
	///	</summary>
	size_t binsize() const {
		return m_map.size();
	}

	///	<summary>
	///		Number of Nodes stored in this object.
	///	</summary>
	size_t size() const {
		size_t total_size = 0;
		typename node_map_type::const_iterator iter = m_map.begin();
		for (; iter != m_map.end(); iter++) {
			total_size += iter->second.size();
		}
		return total_size;
	}

	///	<summary>
	///		Convert a Node to a bin index.
	///	</summary>
	std::array<int,3> NodeToBinIndex(
		const NodeType & node
	) const {
		std::array<int,3> arrBin;
		arrBin[0] = static_cast<int>((node.x + 2.123456789101112) / m_bin_width);
		arrBin[1] = static_cast<int>((node.y + 2.123456789101112) / m_bin_width);
		arrBin[2] = static_cast<int>((node.z + 2.123456789101112) / m_bin_width);
		return arrBin;
	}

	///	<summary>
	///		Insert a given <NodeType, DataType> pair.
	///	</summary>
	void insert(
		const value_type & value
	) {
		std::array<int,3> arrBin = NodeToBinIndex(value.first);

		node_map_type_iterator iter = m_map.find(hasher(arrBin));
		if (iter == m_map.end()) {
			iter = m_map.insert(
				std::pair<bin_index_type, bin_contents_type>(
					hasher(arrBin), bin_contents_type())).first;
		}

		iter->second.push_back(value);
	}

	///	<summary>
	///		Find a given Node.
	///	</summary>
	const_iterator find(
		const NodeType & node
	) const {

		bool fExtended = false;

		// Bins to search over
		std::array<int,3> arrBegin = NodeToBinIndex(node);
		std::array<int,3> arrEnd;

		arrEnd[0] = arrBegin[0] + 1;
		arrEnd[1] = arrBegin[1] + 1;
		arrEnd[2] = arrBegin[2] + 1;

		// Calculate the bin range in the x direction
		double dRemainderXL = node.x - static_cast<double>(arrBegin[0]) * m_bin_width;
		if (fabs(dRemainderXL) < m_tolerance) {
			fExtended = true;
			arrBegin[0]--;
		} else if (fabs(dRemainderXL - 1.0) < m_tolerance) {
			fExtended = true;
			arrEnd[0]++;
		}

		// Calculate the bin range in the y direction
		double dRemainderYL = node.y - static_cast<double>(arrBegin[1]) * m_bin_width;
		if (fabs(dRemainderYL) < m_tolerance) {
			fExtended = true;
			arrBegin[1]--;
		} else if (fabs(dRemainderYL - 1.0) < m_tolerance) {
			fExtended = true;
			arrEnd[1]++;
		}

		// Calculate the bin range in the z direction
		double dRemainderZL = node.z - static_cast<double>(arrBegin[2]) * m_bin_width;
		if (fabs(dRemainderZL) < m_tolerance) {
			fExtended = true;
			arrBegin[2]--;
		} else if (fabs(dRemainderZL - 1.0) < m_tolerance) {
			fExtended = true;
			arrEnd[2]++;
		}

		// Search an extended range of bins
		if (fExtended) {
			return _find_in_bin_range(arrBegin, arrEnd, node);

		} else {
			typename node_map_type::const_iterator iter = m_map.find(hasher(arrBegin));
			if (iter != m_map.end()) {
				size_t ix = _bin_contains(iter->second, node);
				if (ix != ix_not_found) {
					const_iterator iter_out(&m_map);
					iter_out.iter_node_map = iter;
					iter_out.iter_bin_contents = iter->second.begin() + ix;
					return iter_out;
				}
			}
			return end();
		}

	}

private:
	///	<summary>
	///		The internal node map data structure.
	///	</summary>
	node_map_type m_map;

	///	<summary>
	///		The internal floating point tolerance for comparing two Nodes.
	///	</summary>
	const double m_tolerance;

	///	<summary>
	///		The internal m_bin_width.
	///	</summary>
	const double m_bin_width;
};

///////////////////////////////////////////////////////////////////////////////

#endif //_NODEMULTIMAP3D_H_

