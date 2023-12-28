#include "sph_kernel.h"

SPHKernel::SPHKernel(fpreal kernel_radius) : h(kernel_radius), h2(h * h), h3(h2 * h), h4(h2 * h2), h5(h2 * h3) {}
StdKernel::StdKernel(fpreal kernelRadius) : SPHKernel(kernelRadius) {}
SpikyKernel::SpikyKernel(fpreal kernelRadius) : SPHKernel(kernelRadius) {}
auto StdKernel::operator()(fpreal distance) const -> fpreal
{
	if (distance * distance >= h2)
		return 0;

	fpreal x = 1.0 - distance * distance / h2;
	return 315.0 / (64.0 * M_PI * h3) * x * x * x;
}
auto StdKernel::first_derivative(fpreal distance) const -> fpreal
{
	if (distance > h)
		return 0;

	fpreal x = 1.0 - distance * distance / h2;
	return -945.0 / (32.0 * M_PI * h5) * distance * x * x;
}
auto StdKernel::second_derivative(fpreal distance) const -> fpreal
{
	if (distance * distance > h2)
		return 0;

	fpreal x = distance * distance / h2;
	return 945.0 / (32.0 * M_PI * h5) * (1 - x) * (5 * x - 1);
}
auto StdKernel::gradient(const UT_Vector3F &point) const -> UT_Vector3F
{
	fpreal dist = point.length();
	if (dist > 0)
		return gradient(dist, point / dist);
	return UT_Vector3F(0, 0, 0);
}
auto StdKernel::gradient(fpreal distance, const UT_Vector3F &direction) const -> UT_Vector3F
{
	return -first_derivative(distance) * direction;
}

auto StdKernel::cubic_kernel_derivative(fpreal distance) const -> fpreal {
    if (distance > h)
        return 0;
    fpreal x = 1.0 - distance / h;
    return -30.0 / (M_PI * h4) * x * x;
}

auto SpikyKernel::operator()(fpreal distance) const -> fpreal
{
	if (distance > h)
		return 0;

	fpreal x = 1.0 - distance / h;
	return 15.0 / (M_PI * h3) * x * x * x;
}
auto SpikyKernel::first_derivative(fpreal distance) const -> fpreal
{
	if (distance > h)
		return 0;

	fpreal x = 1.0 - distance / h;
	return -45.0 / (M_PI * h4) * x * x;
}
auto SpikyKernel::second_derivative(fpreal distance) const -> fpreal
{
	if (distance > h)
		return 0;

	fpreal x = 1.0 - distance / h;
	return 90.0 / (M_PI * h5) * x;
}
auto SpikyKernel::gradient(const UT_Vector3F &point) const -> UT_Vector3F
{
	fpreal dist = point.length();
	if (dist > 0)
		return gradient(dist, point / dist);
	return UT_Vector3F(0, 0, 0);
}
auto SpikyKernel::gradient(fpreal distance, const UT_Vector3F &direction) const -> UT_Vector3F
{
	return -first_derivative(distance) * direction;
}

