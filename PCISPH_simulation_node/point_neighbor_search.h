#include <SIM/SIM_Utils.h>
#include "size.h"

class PointNeighborSearch3
{
public:
	using ForEachNearbyPointFunc = std::function<void(size_t, const UT_Vector3 &)>;

	virtual void for_each_nearby_point(const UT_Vector3 &origin, const ForEachNearbyPointFunc &callback) = 0;
	virtual auto has_nearby_point(const UT_Vector3 &origin) const -> bool = 0;
	virtual void build(const std::vector<UT_Vector3> &points) = 0;

public:
	explicit PointNeighborSearch3(fpreal radius) : _radius(radius) {}

protected:
	fpreal _radius;
};


class PointSimpleListSearch3 final : public PointNeighborSearch3
{
public:
	void for_each_nearby_point(const UT_Vector3 &origin, const ForEachNearbyPointFunc &callback) final;
	auto has_nearby_point(const UT_Vector3 &origin) const -> bool final;
	void build(const std::vector<UT_Vector3> &points) final;
	explicit PointSimpleListSearch3(fpreal radius) : PointNeighborSearch3(radius) {}

private:
	std::vector<UT_Vector3> _points;
};


class PointHashGridSearch3 final : public PointNeighborSearch3
{
public:
	void for_each_nearby_point(const UT_Vector3 &origin, const ForEachNearbyPointFunc &callback) final;
	auto has_nearby_point(const UT_Vector3 &origin) const -> bool final;
	void build(const std::vector<UT_Vector3> &points) final;
	void add_point(const UT_Vector3 &point);
	explicit PointHashGridSearch3(fpreal radius) : PointNeighborSearch3(radius), _grid_spacing(2 * radius) {}

private:
	auto _get_hash_key_from_position(const UT_Vector3 &position) const -> size_t;
	auto _get_nearby_keys(const UT_Vector3 &position) const -> std::array<size_t, 8>;
	auto _get_hash_key_from_bucket_index(const UT_Vector3i &index) const -> size_t;
	auto _get_bucket_index(const UT_Vector3 &position) const -> UT_Vector3i;

private:
	fpreal _grid_spacing = 0.04; // 2 * kernel radius
	mSize3 _resolution = {10, 10, 10};
	std::vector<UT_Vector3> _points;
	std::vector<std::vector<size_t>> _buckets;
};

class PointParallelHashGridSearch3 final : public PointNeighborSearch3
{
public:
	void for_each_nearby_point(const UT_Vector3 &origin, const ForEachNearbyPointFunc &callback) final;
	auto has_nearby_point(const UT_Vector3 &origin) const -> bool final;
	void build(const std::vector<UT_Vector3> &points) final;
	explicit PointParallelHashGridSearch3(fpreal radius) : PointNeighborSearch3(radius), _grid_spacing(2 * radius) {}

private:
	auto _get_hash_key_from_position(const UT_Vector3 &position) const -> size_t;
	auto _get_nearby_keys(const UT_Vector3 &position) const -> std::array<size_t, 8>;
	auto _get_hash_key_from_bucket_index(const UT_Vector3i &index) const -> size_t;
	auto _get_bucket_index(const UT_Vector3 &position) const -> UT_Vector3i;

private:
	fpreal _grid_spacing = 0.04;
	mSize3 _resolution = {10, 10, 10};
	std::vector<UT_Vector3> _points;
	std::vector<size_t> _keys;
	std::vector<size_t> _start_index_table;
	std::vector<size_t> _end_index_table;
	std::vector<size_t> _sorted_indices;
};

using PointNeighborSearch3Ptr = std::shared_ptr<PointNeighborSearch3>;
using PointSimpleListSearch3Ptr = std::shared_ptr<PointSimpleListSearch3>;
using PointHashGridSearch3Ptr = std::shared_ptr<PointHashGridSearch3>;