#include "point_neighbor_search.h"

// ==================== PointSimpleListSearch3 ====================
void HinaPE::PointSimpleListSearch3::for_each_nearby_point(const UT_Vector3 &origin, const HinaPE::PointNeighborSearch3::ForEachNearbyPointFunc &callback)
{
	fpreal radius_squared = _radius * _radius;

	for (int i = 0; i < _points.size(); ++i)
	{
		UT_Vector3 r = _points[i] - origin;
		fpreal distance_squared = r.dot(r);
		if (distance_squared <= radius_squared)
			callback(i, _points[i]);
	}
}

auto HinaPE::PointSimpleListSearch3::has_nearby_point(const UT_Vector3 &origin) const -> bool
{
	fpreal radius_squared = _radius * _radius;

	return std::any_of(_points.begin(), _points.end(), [radius_squared, origin](const UT_Vector3 &point)
	{
		UT_Vector3 r = point - origin;
		fpreal distance_squared = r.dot(r);
		if (distance_squared <= radius_squared)
			return true;
		return false;
	});
}

void HinaPE::PointSimpleListSearch3::build(const std::vector<UT_Vector3> &points)
{
	_points.resize(points.size());
	std::copy(points.begin(), points.end(), _points.begin());
}


// ==================== PointHashGridSearch3 ====================
void HinaPE::PointHashGridSearch3::for_each_nearby_point(const UT_Vector3 &origin, const HinaPE::PointNeighborSearch3::ForEachNearbyPointFunc &callback)
{
	if (_buckets.empty())
		return;

	auto nearby_keys = _get_nearby_keys(origin);
	const fpreal query_radius_squared = _radius * _radius;

	for (auto key: nearby_keys)
	{
		for (auto i: _buckets[key])
		{
			UT_Vector3 r = _points[i] - origin;
			fpreal distance_squared = r.dot(r);
			if (distance_squared <= query_radius_squared)
				callback(i, _points[i]);
		}
	}
}

auto HinaPE::PointHashGridSearch3::has_nearby_point(const UT_Vector3 &origin) const -> bool
{
	if (_buckets.empty())
		return false;

	auto nearby_keys = _get_nearby_keys(origin);
	const fpreal query_radius_squared = _radius * _radius;

	for (auto key: nearby_keys)
	{
		for (auto i: _buckets[key])
		{
			UT_Vector3 r = _points[i] - origin;
			fpreal distance_squared = r.dot(r);
			if (distance_squared <= query_radius_squared)
				return true;
		}
	}

	return false;
}

void HinaPE::PointHashGridSearch3::build(const std::vector<UT_Vector3> &points)
{
	_buckets.clear();
	_points.clear();

	_buckets.resize(_resolution.x * _resolution.y * _resolution.z);
	_points.resize(points.size());

	if (points.empty()) return;

	for (size_t i = 0; i < points.size(); ++i)
	{
		_points[i] = points[i];
		size_t key = _get_hash_key_from_position(points[i]);
		_buckets[key].push_back(i);
	}
}
void HinaPE::PointHashGridSearch3::add_point(const UT_Vector3 &point)
{
	if (_buckets.empty())
	{
		std::vector<UT_Vector3> points;
		points.push_back(point);
		build(points);
	} else
	{
		size_t i = _points.size();
		_points.push_back(point);
		size_t key = _get_hash_key_from_position(point);
		_buckets[key].push_back(i);
	}
}

auto HinaPE::PointHashGridSearch3::_get_hash_key_from_position(const UT_Vector3 &position) const -> size_t { return _get_hash_key_from_bucket_index(_get_bucket_index(position)); }
auto HinaPE::PointHashGridSearch3::_get_nearby_keys(const UT_Vector3 &position) const -> std::array<size_t, 8>
{
	auto origin_index = _get_bucket_index(position);
	std::array<UT_Vector3i, 8> nearby_bucket_indices{};

	for (int i = 0; i < 8; ++i)
		nearby_bucket_indices[i] = origin_index;

	if ((static_cast<fpreal>(origin_index.x()) + static_cast<fpreal>(0.5)) * _grid_spacing <= position.x())
	{
		nearby_bucket_indices[4].x() += 1;
		nearby_bucket_indices[5].x() += 1;
		nearby_bucket_indices[6].x() += 1;
		nearby_bucket_indices[7].x() += 1;
	} else
	{
		nearby_bucket_indices[4].x() -= 1;
		nearby_bucket_indices[5].x() -= 1;
		nearby_bucket_indices[6].x() -= 1;
		nearby_bucket_indices[7].x() -= 1;
	}

	if ((static_cast<fpreal>(origin_index.y()) + static_cast<fpreal>(0.5)) * _grid_spacing <= position.y())
	{
		nearby_bucket_indices[2].y() += 1;
		nearby_bucket_indices[3].y() += 1;
		nearby_bucket_indices[6].y() += 1;
		nearby_bucket_indices[7].y() += 1;
	} else
	{
		nearby_bucket_indices[2].y() -= 1;
		nearby_bucket_indices[3].y() -= 1;
		nearby_bucket_indices[6].y() -= 1;
		nearby_bucket_indices[7].y() -= 1;
	}

	if ((static_cast<fpreal>(origin_index.z()) + static_cast<fpreal>(0.5)) * _grid_spacing <= position.z())
	{
		nearby_bucket_indices[1].z() += 1;
		nearby_bucket_indices[3].z() += 1;
		nearby_bucket_indices[5].z() += 1;
		nearby_bucket_indices[7].z() += 1;
	} else
	{
		nearby_bucket_indices[1].z() -= 1;
		nearby_bucket_indices[3].z() -= 1;
		nearby_bucket_indices[5].z() -= 1;
		nearby_bucket_indices[7].z() -= 1;
	}

	std::array<size_t, 8> res = {};
	for (int i = 0; i < 8; ++i)
		res[i] = _get_hash_key_from_bucket_index(nearby_bucket_indices[i]);

	return res;
}
auto HinaPE::PointHashGridSearch3::_get_hash_key_from_bucket_index(const UT_Vector3i &index) const -> size_t
{
	UT_Vector3i wrapped_index = index;

	wrapped_index.x() = index.x() % static_cast<int>(_resolution.x);
	wrapped_index.y() = index.y() % static_cast<int>(_resolution.y);
	wrapped_index.z() = index.z() % static_cast<int>(_resolution.z);

	if (wrapped_index.x() < 0) wrapped_index.x() += static_cast<int>(_resolution.x);
	if (wrapped_index.y() < 0) wrapped_index.y() += static_cast<int>(_resolution.y);
	if (wrapped_index.z() < 0) wrapped_index.z() += static_cast<int>(_resolution.z);

	return static_cast<size_t>(wrapped_index.z() * _resolution.y + wrapped_index.y()) * _resolution.x + wrapped_index.x();
}
auto HinaPE::PointHashGridSearch3::_get_bucket_index(const UT_Vector3 &position) const -> UT_Vector3i
{
	UT_Vector3i res;
	res.x() = static_cast<int>(std::floor(position.x() / _grid_spacing));
	res.y() = static_cast<int>(std::floor(position.y() / _grid_spacing));
	res.z() = static_cast<int>(std::floor(position.z() / _grid_spacing));
	return res;
}


// ==================== PointParallelHashGridSearch3 ====================
void HinaPE::PointParallelHashGridSearch3::for_each_nearby_point(const UT_Vector3 &origin, const HinaPE::PointNeighborSearch3::ForEachNearbyPointFunc &callback)
{
	auto nearby_keys = _get_nearby_keys(origin);
	const auto query_radius_squared = _radius * _radius;

	for (auto nearby_key: nearby_keys)
	{
		size_t start = _start_index_table[nearby_key];
		size_t end = _end_index_table[nearby_key];
//		|| end == +std::numeric_limits<size_t>::max()
		if (start == +std::numeric_limits<size_t>::max())
			continue;

		for (size_t i = start; i < end; ++i)
		{
			UT_Vector3 direction = _points[i] - origin;
			fpreal distance_squared = direction.length2();
			if (distance_squared <= query_radius_squared)
				callback(_sorted_indices[i], _points[i]);
		}
	}
}
auto HinaPE::PointParallelHashGridSearch3::has_nearby_point(const UT_Vector3 &origin) const -> bool
{
	auto nearby_keys = _get_nearby_keys(origin);
	const fpreal query_radius_squared = _radius * _radius;


	for (auto nearby_key: nearby_keys)
	{
		size_t start = _start_index_table[nearby_key];
		size_t end = _end_index_table[nearby_key];

		if (start == +std::numeric_limits<size_t>::max() || end == +std::numeric_limits<size_t>::max())
			continue;

		for (size_t i = start; i < end; ++i)
		{
			UT_Vector3 direction = _points[i] - origin;
			fpreal distance_squared = direction.length2();
			if (distance_squared <= query_radius_squared)
				return true;
		}
	}

	return false;
}
void HinaPE::PointParallelHashGridSearch3::build(const std::vector<UT_Vector3> &points)
{
	if (points.empty())
		return;

	_points.clear();
	_keys.clear();
	_start_index_table.clear();
	_end_index_table.clear();
	_sorted_indices.clear();

	auto num_points = points.size();
	std::vector<size_t> temp_keys(num_points);
	_start_index_table.resize(_resolution.x * _resolution.y * _resolution.z, +std::numeric_limits<size_t>::max());
	_end_index_table.resize(_resolution.x * _resolution.y * _resolution.z, +std::numeric_limits<size_t>::max());
	_keys.resize(num_points);
	_sorted_indices.resize(num_points);
	_points.resize(num_points);

	// Initialize indices array and generate hash key for each point

    for(int i = 0; i < num_points; ++i)
    {
        _sorted_indices[i] = i;
        _points[i] = points[i];
        temp_keys[i] = _get_hash_key_from_position(points[i]);
    }

	// Sort indices based on hash key
    std::sort(_sorted_indices.begin(), _sorted_indices.end(), [&](size_t a, size_t b)
    {
        return temp_keys[a] < temp_keys[b];
    });

	// Re-order point and key arrays
    for(size_t i = 0; i < num_points; ++i)
    {
        _points[i] = points[_sorted_indices[i]];
        _keys[i] = temp_keys[_sorted_indices[i]];
    }

	// Now _points and _keys are sorted by points' hash key values.
	// Let's fill in start/end index table with _keys.

	// Assume that _keys array looks like:
	// [5|8|8|10|10|10]
	// Then _startIndexTable and _endIndexTable should be like:
	// [.....|0|...|1|..|3|..]
	// [.....|1|...|3|..|6|..]
	//       ^5    ^8   ^10
	// So that _endIndexTable[i] - _startIndexTable[i] is the number points
	// in i-th table bucket.
	_start_index_table[_keys[0]] = 0;
	_end_index_table[_keys[num_points - 1]] = num_points;

    for(size_t i = 1; i < num_points; ++i)
    {
        if (_keys[i] > _keys[i - 1])
        {
            _start_index_table[_keys[i]] = i;
            _end_index_table[_keys[i - 1]] = i;
        }
    }

	size_t sum_num_of_points_per_bucket = 0;
	size_t max_num_of_points_per_bucket = 0;
	size_t num_of_non_empty_buckets = 0;
	for (size_t i = 0; i < _start_index_table.size(); ++i)
	{
		if (_start_index_table[i] != +std::numeric_limits<size_t>::max())
		{
			size_t num_of_points_in_bucket = _end_index_table[i] - _start_index_table[i];
			sum_num_of_points_per_bucket += num_of_points_in_bucket;
			max_num_of_points_per_bucket = std::max(max_num_of_points_per_bucket, num_of_points_in_bucket);
			++num_of_non_empty_buckets;
		}
	}
}
auto HinaPE::PointParallelHashGridSearch3::_get_hash_key_from_position(const UT_Vector3 &position) const -> size_t { return _get_hash_key_from_bucket_index(_get_bucket_index(position)); }
auto HinaPE::PointParallelHashGridSearch3::_get_nearby_keys(const UT_Vector3 &position) const -> std::array<size_t, 8>
{
	auto origin_index = _get_bucket_index(position);
	std::array<UT_Vector3i, 8> nearby_bucket_indices{};

	for (int i = 0; i < 8; ++i)
		nearby_bucket_indices[i] = origin_index;

	if ((static_cast<fpreal>(origin_index.x()) + static_cast<fpreal>(0.5)) * _grid_spacing <= position.x())
	{
		nearby_bucket_indices[4].x() += 1;
		nearby_bucket_indices[5].x() += 1;
		nearby_bucket_indices[6].x() += 1;
		nearby_bucket_indices[7].x() += 1;
	} else
	{
		nearby_bucket_indices[4].x() -= 1;
		nearby_bucket_indices[5].x() -= 1;
		nearby_bucket_indices[6].x() -= 1;
		nearby_bucket_indices[7].x() -= 1;
	}

	if ((static_cast<fpreal>(origin_index.y()) + static_cast<fpreal>(0.5)) * _grid_spacing <= position.y())
	{
		nearby_bucket_indices[2].y() += 1;
		nearby_bucket_indices[3].y() += 1;
		nearby_bucket_indices[6].y() += 1;
		nearby_bucket_indices[7].y() += 1;
	} else
	{
		nearby_bucket_indices[2].y() -= 1;
		nearby_bucket_indices[3].y() -= 1;
		nearby_bucket_indices[6].y() -= 1;
		nearby_bucket_indices[7].y() -= 1;
	}

	if ((static_cast<fpreal>(origin_index.z()) + static_cast<fpreal>(0.5)) * _grid_spacing <= position.z())
	{
		nearby_bucket_indices[1].z() += 1;
		nearby_bucket_indices[3].z() += 1;
		nearby_bucket_indices[5].z() += 1;
		nearby_bucket_indices[7].z() += 1;
	} else
	{
		nearby_bucket_indices[1].z() -= 1;
		nearby_bucket_indices[3].z() -= 1;
		nearby_bucket_indices[5].z() -= 1;
		nearby_bucket_indices[7].z() -= 1;
	}

	std::array<size_t, 8> res = {};
	for (int i = 0; i < 8; ++i)
		res[i] = _get_hash_key_from_bucket_index(nearby_bucket_indices[i]);

	return res;
}
auto HinaPE::PointParallelHashGridSearch3::_get_hash_key_from_bucket_index(const UT_Vector3i &index) const -> size_t
{
	UT_Vector3i wrapped_index = index;

	wrapped_index.x() = index.x() % static_cast<int>(_resolution.x);
	wrapped_index.y() = index.y() % static_cast<int>(_resolution.y);
	wrapped_index.z() = index.z() % static_cast<int>(_resolution.z);

	if (wrapped_index.x() < 0) wrapped_index.x() += static_cast<int>(_resolution.x);
	if (wrapped_index.y() < 0) wrapped_index.y() += static_cast<int>(_resolution.y);
	if (wrapped_index.z() < 0) wrapped_index.z() += static_cast<int>(_resolution.z);

	return static_cast<size_t>(wrapped_index.z() * _resolution.y + wrapped_index.y()) * _resolution.x + wrapped_index.x();
}
auto HinaPE::PointParallelHashGridSearch3::_get_bucket_index(const UT_Vector3 &position) const -> UT_Vector3i
{
	UT_Vector3i res;
	res.x() = static_cast<int>(std::floor(position.x() / _grid_spacing));
	res.y() = static_cast<int>(std::floor(position.y() / _grid_spacing));
	res.z() = static_cast<int>(std::floor(position.z() / _grid_spacing));
	return res;
}
