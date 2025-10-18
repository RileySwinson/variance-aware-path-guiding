/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/render/testcase.h>
#include <mitsuba/core/quad.h>
#include <boost/bind.hpp>

MTS_NAMESPACE_BEGIN

class TestQuadrature : public TestCase {
public:
	MTS_BEGIN_TESTCASE()
	MTS_DECLARE_TEST(test01_quad)
	MTS_DECLARE_TEST(test02_nD_01)
	MTS_DECLARE_TEST(test03_nD_02)
	MTS_DECLARE_TEST(test04_gaussLegendre_even)
	MTS_DECLARE_TEST(test05_gaussLegendre_odd)
	MTS_DECLARE_TEST(test06_gaussLobatto_even)
	MTS_DECLARE_TEST(test07_gaussLobatto_odd)
	MTS_END_TESTCASE()

	float testF(float t) const {
		return std::sin(t);
	}

	void testF2(const float *in, float *out) const {
		*out = std::sin(*in);
	}

	inline float gauss3(Vector x, float stddev) const {
		return math::fastexp(-0.5f * dot(x, x)/stddev)/(std::pow(2*M_PI * stddev, (float) 3 / (float) 2));
	}

	void testF3(size_t nPoints, const float *in, float *out) const {
		for (size_t i=0; i<nPoints; ++i) {
			Vector v(in[0], in[1], in[2]);
			float weight =
				(1 + v.x*v.x) / std::pow(1-v.x*v.x, 2) *
				(1 + v.y*v.y) / std::pow(1-v.y*v.y, 2) *
				(1 + v.z*v.z) / std::pow(1-v.z*v.z, 2);
			v = Vector(v.x / (1-v.x*v.x), v.y / (1-v.y*v.y),
					   v.z / (1-v.z*v.z));
			out[0*nPoints + i] = weight * gauss3(v, 0.1);
			out[1*nPoints + i] = weight * gauss3(v, 1);
			in += 3;
		}
	}

	void test01_quad() {
		GaussLobattoIntegrator quad(1024, 0, 1e-5f);
		size_t evals;
		float result = quad.integrate(boost::bind(
			&TestQuadrature::testF, this, _1), 0, 10, &evals);
		float ref = 2 * std::pow(std::sin((float) 5.0f), (float) 2.0f);
		Log(EInfo, "test01_quad(): used " SIZE_T_FMT " function evaluations", evals);
		assertEqualsEpsilon(result, ref, 1e-5f);
	}

	void test02_nD_01() {
		NDIntegrator quad(1, 1, 1024, 0, 1e-5f);
		float min = 0, max = 10, result, err;
		size_t evals;
		assertTrue(quad.integrate(boost::bind(
			&TestQuadrature::testF2, this, _1, _2), &min, &max, &result, &err, &evals) == NDIntegrator::ESuccess);
		float ref = 2 * std::pow(std::sin(5.0f), 2.0f);
		Log(EInfo, "test02_nD_01(): used " SIZE_T_FMT " function evaluations, "
				"error=%f", evals, err);
		assertEqualsEpsilon(result, ref, 1e-5f);
	}

	void test03_nD_02() {
		NDIntegrator quad(2, 3, 1000000, 0, 1e-5f);
		size_t evals;
		float min[3] = { -1, -1, -1 } , max[3] = { 1, 1, 1 }, result[2], err[2];
		assertTrue(quad.integrateVectorized(boost::bind(
			&TestQuadrature::testF3, this, _1, _2, _3), min, max, result, err, &evals) == NDIntegrator::ESuccess);
		Log(EInfo, "test02_nD_02(): used " SIZE_T_FMT " function evaluations, "
				"error=[%f, %f]", evals, err[0], err[1]);
		assertEqualsEpsilon(result[0], 1.0f, 1e-5f);
		assertEqualsEpsilon(result[1], 1.0f, 1e-5f);
	}

	void test04_gaussLegendre_even() {
		float nodes[4], weights[4];
		gaussLegendre(4, nodes, weights);

		assertEqualsEpsilon(nodes[0], (float) (-1/35.0 * std::sqrt(525+70*std::sqrt(30.0))), 1e-8f);
		assertEqualsEpsilon(nodes[1], (float) (-1/35.0 * std::sqrt(525-70*std::sqrt(30.0))), 1e-8f);
		assertEqualsEpsilon(nodes[2], (float) ( 1/35.0 * std::sqrt(525-70*std::sqrt(30.0))), 1e-8f);
		assertEqualsEpsilon(nodes[3], (float) ( 1/35.0 * std::sqrt(525+70*std::sqrt(30.0))), 1e-8f);

		assertEqualsEpsilon(weights[0], (float) (1.0/36.0 * (18-std::sqrt(30.0))), 1e-8f);
		assertEqualsEpsilon(weights[1], (float) (1.0/36.0 * (18+std::sqrt(30.0))), 1e-8f);
		assertEqualsEpsilon(weights[2], (float) (1.0/36.0 * (18+std::sqrt(30.0))), 1e-8f);
		assertEqualsEpsilon(weights[3], (float) (1.0/36.0 * (18-std::sqrt(30.0))), 1e-8f);
	}

	void test05_gaussLegendre_odd() {
		float nodes[5], weights[5];
		gaussLegendre(5, nodes, weights);

		assertEqualsEpsilon(nodes[0], (float) (-1/21.0 * std::sqrt(245+14*std::sqrt(70.0))), 1e-8f);
		assertEqualsEpsilon(nodes[1], (float) (-1/21.0 * std::sqrt(245-14*std::sqrt(70.0))), 1e-8f);
		assertEqualsEpsilon(nodes[2], 0.0, 1e-8f);
		assertEqualsEpsilon(nodes[3], (float) ( 1/21.0 * std::sqrt(245-14*std::sqrt(70.0))), 1e-8f);
		assertEqualsEpsilon(nodes[4], (float) ( 1/21.0 * std::sqrt(245+14*std::sqrt(70.0))), 1e-8f);

		assertEqualsEpsilon(weights[0], (float) (1.0/900.0 * (322-13*std::sqrt(70.0))), 1e-8f);
		assertEqualsEpsilon(weights[1], (float) (1.0/900.0 * (322+13*std::sqrt(70.0))), 1e-8f);
		assertEqualsEpsilon(weights[2], (float) 128.0/225.0, 1e-8f);
		assertEqualsEpsilon(weights[3], (float) (1.0/900.0 * (322+13*std::sqrt(70.0))), 1e-8f);
		assertEqualsEpsilon(weights[4], (float) (1.0/900.0 * (322-13*std::sqrt(70.0))), 1e-8f);
	}

	void test06_gaussLobatto_even() {
		float nodes[6], weights[6];
		gaussLobatto(6, nodes, weights);

		assertEqualsEpsilon(nodes[0], -1.0f, 1e-8f);
		assertEqualsEpsilon(nodes[1], (float) -std::sqrt(1/21.0 * (7+2*std::sqrt(7.0))), 1e-8f);
		assertEqualsEpsilon(nodes[2], (float) -std::sqrt(1/21.0 * (7-2*std::sqrt(7.0))), 1e-8f);
		assertEqualsEpsilon(nodes[3], (float) std::sqrt(1/21.0 * (7-2*std::sqrt(7.0))), 1e-8f);
		assertEqualsEpsilon(nodes[4], (float) std::sqrt(1/21.0 * (7+2*std::sqrt(7.0))), 1e-8f);
		assertEqualsEpsilon(nodes[5], 1.0f, 1e-8f);

		assertEqualsEpsilon(weights[0], (float) (1.0/15.0), 1e-8f);
		assertEqualsEpsilon(weights[1], (float) ((14 - std::sqrt(7.0))/30.0), 1e-8f);
		assertEqualsEpsilon(weights[2], (float) ((14 + std::sqrt(7.0))/30.0), 1e-8f);
		assertEqualsEpsilon(weights[3], (float) ((14 + std::sqrt(7.0))/30.0), 1e-8f);
		assertEqualsEpsilon(weights[4], (float) ((14 - std::sqrt(7.0))/30.0), 1e-8f);
		assertEqualsEpsilon(weights[5], (float) (1.0/15.0), 1e-8f);
	}

	void test07_gaussLobatto_odd() {
		float nodes[5], weights[5];
		gaussLobatto(5, nodes, weights);

		assertEqualsEpsilon(nodes[0], -1.0f, 1e-8f);
		assertEqualsEpsilon(nodes[1], (float) -(std::sqrt(21.0)/7.0), 1e-8f);
		assertEqualsEpsilon(nodes[2], 0.0f, 1e-8f);
		assertEqualsEpsilon(nodes[3], (float) (std::sqrt(21.0)/7.0), 1e-8f);
		assertEqualsEpsilon(nodes[4], 1.0f, 1e-8f);

		assertEqualsEpsilon(weights[0], (float) (1.0/10.0), 1e-8f);
		assertEqualsEpsilon(weights[1], (float) (49.0/90.0), 1e-8f);
		assertEqualsEpsilon(weights[2], (float) (32.0/45.0), 1e-8f);
		assertEqualsEpsilon(weights[3], (float) (49.0/90.0), 1e-8f);
		assertEqualsEpsilon(weights[4], (float) (1.0/10.0), 1e-8f);
	}
};

MTS_EXPORT_TESTCASE(TestQuadrature, "Testcase for quadrature routines")
MTS_NAMESPACE_END
