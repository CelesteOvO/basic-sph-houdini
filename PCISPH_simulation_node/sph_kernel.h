#include <SIM/SIM_Utils.h>

struct SPHKernel
{
public:
	virtual auto operator()(fpreal distance) const -> fpreal = 0;
	virtual auto first_derivative(fpreal distance) const -> fpreal = 0;
	virtual auto second_derivative(fpreal distance) const -> fpreal = 0;
	virtual auto gradient(const UT_Vector3F &point) const -> UT_Vector3F = 0;
	virtual auto gradient(fpreal distance, const UT_Vector3F &direction) const -> UT_Vector3F = 0;
	explicit SPHKernel(fpreal kernel_radius);

public:
	fpreal h, h2, h3, h4, h5;
};

struct StdKernel : public SPHKernel
{
public:
	auto operator()(fpreal distance) const -> fpreal final;
	auto first_derivative(fpreal distance) const -> fpreal final;
	auto second_derivative(fpreal distance) const -> fpreal final;
	auto gradient(const UT_Vector3F &point) const -> UT_Vector3F final;
	auto gradient(fpreal distance, const UT_Vector3F &direction) const -> UT_Vector3F final;
    auto cubic_kernel_derivative(fpreal distance) const -> fpreal ; // temp here
	explicit StdKernel(fpreal kernelRadius);
};

struct SpikyKernel : public SPHKernel
{
public:
	auto operator()(fpreal distance) const -> fpreal final;
	auto first_derivative(fpreal distance) const -> fpreal final;
	auto second_derivative(fpreal distance) const -> fpreal final;
	auto gradient(const UT_Vector3F &point) const -> UT_Vector3F final;
	auto gradient(fpreal distance, const UT_Vector3F &direction) const -> UT_Vector3F final;
	explicit SpikyKernel(fpreal kernelRadius);
};

