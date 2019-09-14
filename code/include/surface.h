#pragma once

#include <vector.h>
#include <random.h>
#include <cfloat>

#define INV_M_PI		0.31830988618379067153f /* 1/pi */
#define SQRT_M_PI		1.77245385090551602729f /* sqrt(pi) */
#define SQRT_2			1.41421356237309504880f /* sqrt(2) */
#define INV_SQRT_M_PI	0.56418958354775628694f /* 1/sqrt(pi) */
#define INV_2_SQRT_M_PI	0.28209479177387814347f /* 0.5/sqrt(pi) */
#define INV_SQRT_2_M_PI 0.3989422804014326779f /* 1/sqrt(2*pi) */
#define INV_SQRT_2		0.7071067811865475244f /* 1/sqrt(2) */

const float sign( const float a )
{
    if( a > 0.0f )
    {
        return 1.0f;
    }
    else
    {
        return -1.0f;
    }
}

static bool IsFiniteNumber(float x) 
{
	return (x <= FLT_MAX && x >= -FLT_MAX); 
} 

static float erfinv(float x)
{
float w, p;
w = - logf((1.0f-x)*(1.0f+x));
if ( w < 5.000000f ) {
w = w - 2.500000f;
p = 2.81022636e-08f;
p = 3.43273939e-07f + p*w;
p = -3.5233877e-06f + p*w;
p = -4.39150654e-06f + p*w;
p = 0.00021858087f + p*w;
p = -0.00125372503f + p*w;
p = -0.00417768164f + p*w;
p = 0.246640727f + p*w;
p = 1.50140941f + p*w;
}
else {
w = sqrtf(w) - 3.000000f;
p = -0.000200214257f;
p = 0.000100950558f + p*w;
p = 0.00134934322f + p*w;
p = -0.00367342844f + p*w;
p = 0.00573950773f + p*w;
p = -0.0076224613f + p*w;
p = 0.00943887047f + p*w;
p = 1.00167406f + p*w;
p = 2.83297682f + p*w;
}
return p*x;
}

//////////////////////////////////////////////////////////////////////////////////
// rough dielectric surface sampling [Heitz et al. 2016]
//  "Multiple-Scattering Microfacet BSDFs with the Smith Model"
//////////////////////////////////////////////////////////////////////////////////

class MicrosurfaceSlope
{
public:
	MicrosurfaceSlope(const float alpha_x=1.0f, const float alpha_y=1.0f)
		: m_alpha_x(alpha_x), m_alpha_y(alpha_y)
	{}

public:
	// roughness
	const float m_alpha_x, m_alpha_y;
	// projected roughness in wi
	float alpha_i(const Vector3& wi) const; 

public:
	// distribution of normals (NDF)	
	float D(const Vector3& wm) const; 
	// distribution of visible normals (VNDF)
	float D_wi(const Vector3& wi, const Vector3& wm) const; 
	// sample the VNDF
	Vector3 sampleD_wi(const Vector3& wi, const float U1, const float U2) const;

public:
	// distribution of slopes
	virtual float P22(const float slope_x, const float slope_y) const=0; 
	// Smith's Lambda function
	virtual float Lambda(const Vector3& wi) const=0;
	// projected area towards incident direction
	virtual float projectedArea(const Vector3& wi) const=0;
	// sample the distribution of visible slopes with alpha=1.0
	virtual Vector2 sampleP22_11(const float theta_i, const float U1, const float U2) const=0;
};

float MicrosurfaceSlope::D(const Vector3& wm) const {
	if( wm.z <= 0.0f)
		return 0.0f;

	// slope of wm
	const float slope_x = -wm.x/wm.z;
	const float slope_y = -wm.y/wm.z;

	// value
	const float value = P22(slope_x, slope_y) / (wm.z*wm.z*wm.z*wm.z);
	return value;
}

Vector3 MicrosurfaceSlope::sampleD_wi(const Vector3& wi, const float U1, const float U2) const {

	float l_alpha_x( m_alpha_x );
	float l_alpha_y( m_alpha_y );

	// stretch to match configuration with alpha=1.0	
	const Vector3 wi_11 = Normalize(Vector3(l_alpha_x * wi.x, l_alpha_y * wi.y, wi.z));

	// sample visible slope with alpha=1.0
	Vector2 slope_11 = sampleP22_11(acosf(wi_11.z), U1, U2);

	// align with view direction
	const float phi = atan2(wi_11.y, wi_11.x);
	Vector2 slope(cosf(phi)*slope_11.x - sinf(phi)*slope_11.y, sinf(phi)*slope_11.x + cos(phi)*slope_11.y);

	// stretch back
	slope.x *= l_alpha_x;
	slope.y *= l_alpha_y;

	// if numerical instability
	if( (slope.x != slope.x) || !IsFiniteNumber(slope.x) ) 
	{
		if(wi.z > 0) return Vector3(0.0f,0.0f,1.0f);
		else return Normalize(Vector3(wi.x, wi.y, 0.0f));
	}

	// compute normal
	const Vector3 wm = Normalize(Vector3(-slope.x, -slope.y, 1.0f));
	return wm;
}

float MicrosurfaceSlope::alpha_i(const Vector3& wi) const
{
	const float invSinTheta2 = 1.0f / (1.0f - wi.z*wi.z);
	const float cosPhi2 = wi.x*wi.x*invSinTheta2;
	const float sinPhi2 = wi.y*wi.y*invSinTheta2;
	const float alpha_i = sqrtf( cosPhi2*m_alpha_x*m_alpha_x + sinPhi2*m_alpha_y*m_alpha_y );
	return alpha_i;
}

float MicrosurfaceSlope::D_wi(const Vector3& wi, const Vector3& wm) const {
	if( wm.z <= 0.0f)
		return 0.0f;

	// normalization coefficient
	const float projectedarea = projectedArea(wi);
	if(projectedarea == 0)
		return 0;
	const float c = 1.0f / projectedarea;

	// value
	const float value = c * std::max(0.0, Dot(wi, wm)) * D(wm);
	return value;
}

class MicrosurfaceSlopeBeckmann : public MicrosurfaceSlope
{
public:
	MicrosurfaceSlopeBeckmann(const float alpha_x=1.0f, const float alpha_y=1.0f)
		: MicrosurfaceSlope(alpha_x, alpha_y)
	{}

	// distribution of slopes
	virtual float P22(const float slope_x, const float slope_y) const; 
	// Smith's Lambda function
	virtual float Lambda(const Vector3& wi) const;
	// projected area towards incident direction
	virtual float projectedArea(const Vector3& wi) const;
	// sample the distribution of visible slopes with alpha=1.0
	virtual Vector2 sampleP22_11(const float theta_i, const float U1, const float U2) const;
};

float MicrosurfaceSlopeBeckmann::P22(const float slope_x, const float slope_y) const
{
	const float value = 1.0f / (M_PI * m_alpha_x * m_alpha_y) * expf(-slope_x*slope_x/(m_alpha_x*m_alpha_x) - slope_y*slope_y/(m_alpha_y*m_alpha_y) );
	return value;
}

float MicrosurfaceSlopeBeckmann::projectedArea(const Vector3& wi) const
{
	if(wi.z > 0.9999f)
		return 1.0f;
	if(wi.z < -0.9999f)
		return 0.0f;

	// a
	const float alphai = alpha_i(wi);
	const float theta_i = acosf(wi.z);
	const float sin_theta_i = sinf(theta_i);
	const float a = 1.0f/tanf(theta_i)/alphai;

	// value
	const float value = 0.5f*((float)erf(a) + 1.0f)*wi.z + INV_2_SQRT_M_PI * alphai * sinf(theta_i) * expf(-a*a);

	return value;
}

float MicrosurfaceSlopeBeckmann::Lambda(const Vector3& wi) const
{
	if(wi.z > 0.9999f)
		return 0.0f;
	if(wi.z < -0.9999f)
		return -1.0f;

	// a
	const float theta_i = acosf(wi.z);
	const float a = 1.0f/tanf(theta_i)/alpha_i(wi);

	// value
	const float value = 0.5f*((float)erf(a) - 1.0f) + INV_2_SQRT_M_PI / a * expf(-a*a);

	return value;
}

Vector2 MicrosurfaceSlopeBeckmann::sampleP22_11(const float theta_i, const float U, const float U_2) const
{
	Vector2 slope;

	if(theta_i < 0.0001f)
	{
		const float r = sqrtf(-logf(U));
		const float phi = 6.28318530718f * U_2;
		slope.x = r * cosf(phi);
		slope.y = r * sinf(phi);
		return slope;
	}

	// constant
	const float sin_theta_i = sinf(theta_i);
	const float cos_theta_i = cosf(theta_i);

	// slope associated to theta_i
	const float slope_i = cos_theta_i/sin_theta_i;

	// projected area
	const float a = cos_theta_i/sin_theta_i;	
	const float projectedarea = 0.5f*((float)erf(a) + 1.0f)*cos_theta_i + INV_2_SQRT_M_PI * sin_theta_i * expf(-a*a);
	if(projectedarea < 0.0001f || projectedarea!=projectedarea)
		return Vector2(0,0);
	// VNDF normalization factor
	const float c = 1.0f / projectedarea;

	// search 
	float erf_min = -0.9999f;
	float erf_max = std::max(erf_min, (float)erf(slope_i));
	float erf_current = 0.5f * (erf_min+erf_max);

	while(erf_max-erf_min > 0.00001f)
	{
		if (!(erf_current >= erf_min && erf_current <= erf_max))
			erf_current = 0.5f * (erf_min + erf_max);

		// evaluate slope
		const float slope = erfinv(erf_current);

		// CDF
		const float CDF = (slope>=slope_i) ? 1.0f : c * (INV_2_SQRT_M_PI*sin_theta_i*expf(-slope*slope) + cos_theta_i*(0.5f+0.5f*(float)erf(slope)));
		const float diff = CDF - U;

		// test estimate
		if( abs(diff) < 0.00001f )
			break;

		// update bounds
		if(diff > 0.0f)
		{
			if(erf_max == erf_current)
				break;
			erf_max = erf_current;
		}
		else
		{
			if(erf_min == erf_current)
				break;
			erf_min = erf_current;
		}

		// update estimate
		const float derivative = 0.5f*c*cos_theta_i - 0.5f*c*sin_theta_i * slope;
		erf_current -= diff/derivative;
	}

	slope.x = erfinv(std::min(erf_max, std::max(erf_min, erf_current)));
	slope.y = erfinv(2.0f*U_2-1.0f);
	return slope;
}

class MicrosurfaceSlopeGGX : public MicrosurfaceSlope
{
public:
	MicrosurfaceSlopeGGX(const float alpha_x=1.0f, const float alpha_y=1.0f)
		: MicrosurfaceSlope(alpha_x, alpha_y)
	{}

	// distribution of slopes
	virtual float P22(const float slope_x, const float slope_y) const; 
	// Smith's Lambda function
	virtual float Lambda(const Vector3& wi) const;
	// projected area towards incident direction
	virtual float projectedArea(const Vector3& wi) const;
	// sample the distribution of visible slopes with alpha=1.0
	virtual Vector2 sampleP22_11(const float theta_i, const float U1, const float U2) const;
};

float MicrosurfaceSlopeGGX::P22(const float slope_x, const float slope_y) const
{
	const float tmp = 1.0f + slope_x*slope_x/(m_alpha_x*m_alpha_x) + slope_y*slope_y/(m_alpha_y*m_alpha_y);
	const float value = 1.0f / (M_PI * m_alpha_x * m_alpha_y) / (tmp * tmp);
	return value;
}

float MicrosurfaceSlopeGGX::Lambda(const Vector3& wi) const
{
	if(wi.z > 0.9999f)
		return 0.0f;
	if(wi.z < -0.9999f)
		return -1.0f;

	// a
	const float theta_i = acosf(wi.z);
	const float a = 1.0f/tanf(theta_i)/alpha_i(wi);

	// value
	const float value = 0.5f*(-1.0f + sign(a) * sqrtf(1 + 1/(a*a)));

	return value;
}

Vector2 MicrosurfaceSlopeGGX::sampleP22_11(const float theta_i, const float U, const float U_2) const
{
	Vector2 slope;

	if(theta_i < 0.0001f)
	{
		const float r = sqrtf(U/(1.0f-U));
		const float phi = 6.28318530718f * U_2;
		slope.x = r * cosf(phi);
		slope.y = r * sinf(phi);
		return slope;
	}

	// constant
	const float sin_theta_i = sinf(theta_i);
	const float cos_theta_i = cosf(theta_i);
	const float tan_theta_i = sin_theta_i/cos_theta_i;

	// slope associated to theta_i
	const float slope_i = cos_theta_i/sin_theta_i;

	// projected area
	const float projectedarea = 0.5f * (cos_theta_i + 1.0f);
	if(projectedarea < 0.0001f || projectedarea!=projectedarea)
		return Vector2(0,0);
	// normalization coefficient
	const float c = 1.0f / projectedarea;

	const float A = 2.0f*U/cos_theta_i/c - 1.0f;
	const float B = tan_theta_i;
	const float tmp = 1.0f / (A*A-1.0f);

	const float D = sqrtf(std::max(0.0f, B*B*tmp*tmp - (A*A-B*B)*tmp));
	const float slope_x_1 = B*tmp - D;
	const float slope_x_2 = B*tmp + D;
	slope.x = (A < 0.0f || slope_x_2 > 1.0f/tan_theta_i) ? slope_x_1 : slope_x_2;
	slope.y = Sqrt(-1 - slope.x*slope.x + (1 + slope.x*slope.x)/Power(1 - U_2, 0.6666666666666666))*sin(2*Pi*RandomReal());

	return slope;
}

float MicrosurfaceSlopeGGX::projectedArea(const Vector3& wi) const
{
	if(wi.z > 0.9999f)
		return 1.0f;
	if( wi.z < -0.9999f)
		return 0.0f;

	// a
	const float theta_i = acosf(wi.z);
	const float sin_theta_i = sinf(theta_i);

	const float alphai = alpha_i(wi);

	// value
	const float value = 0.5f * (wi.z + sqrtf(wi.z*wi.z + sin_theta_i*sin_theta_i*alphai*alphai));

	return value;
}

class Microsurface 
{
public:
	// slope distribution
	const MicrosurfaceSlope* m_microsurfaceslope; 

public:

	Microsurface(const bool slope_beckmann, // Beckmann or GGX slope distribution
				const float alpha_x,
				const float alpha_y ) : 
		m_microsurfaceslope((slope_beckmann) ? 
		  static_cast<MicrosurfaceSlope*>(new MicrosurfaceSlopeBeckmann(alpha_x, alpha_y)) 
		: static_cast<MicrosurfaceSlope*>(new MicrosurfaceSlopeGGX(alpha_x, alpha_y)))
	{}
		
	~Microsurface() 
	{
		
	}

	// evaluate BSDF with a random walk (stochastic but unbiased)
	// scatteringOrder=0 --> contribution from all scattering events
	// scatteringOrder=1 --> contribution from 1st bounce only
	// scatteringOrder=2 --> contribution from 2nd bounce only, etc..
	//virtual float eval(const Vector3& wi, const Vector3& wo, const int scatteringOrder=0) const; 

	// evaluate BSDF with bidirectional algorithm for variance reduction
	//virtual float evalBidir(const Vector3& wi, const Vector3& wo, const int scatteringOrder=0) const;
	//virtual float evalMIS(const Vector3& wi, const Vector3& wo, const int scatteringOrder=0) const; 

	// sample BSDF with a random walk
	// scatteringOrder is set to the number of bounces computed for this sample
	virtual Vector3 sample(const Vector3& wi, int& scatteringOrder) const;
	Vector3 sample(const Vector3& wi) const {int scatteringOrder; return sample(wi, scatteringOrder);}

public:
	// masking function
	float G_1(const Vector3& wi) const;
	// masking function at height h0
	float G_1(const Vector3& wi, const float h0) const;
	// sample height in outgoing direction
	float sampleHeight(const Vector3& wo, const float h0, const float U) const;

public:
	// evaluate local phase function 
	//virtual float evalPhaseFunction(const Vector3& wi, const Vector3& wo) const=0;
	// sample local phase function
	virtual Vector3 samplePhaseFunction(const Vector3& wi) const=0; 

	// evaluate BSDF limited to single scattering 
	// this is in average equivalent to eval(wi, wo, 1);
	//virtual float evalSingleScattering(const Vector3& wi, const Vector3& wo) const=0; 
};

Vector3 Microsurface::sample(const Vector3& wi, int& scatteringOrder) const
{
	// init
	Vector3 wr = -wi;
	float hr = 0.0f;
	
	// random walk
	scatteringOrder = 0;	
	while(true)
	{
		// next height
		float U = RandomReal();
		hr = sampleHeight(wr, hr, U);

		// leave the microsurface?
		if( hr >= 0.0f )
			break;
		else
			scatteringOrder++;

		// next direction
		wr = samplePhaseFunction(-wr);

		// if NaN (should not happen, just in case)
		if( (hr != hr) || (wr.z != wr.z) ) 
			return Vector3(0,0,1);
	}

	return wr;
}

float Microsurface::G_1(const Vector3& wi) const
{
	if(wi.z > 0.9999f)
		return 1.0f;
	if(wi.z <= 0.0f)
		return 0.0f;

	// Lambda
	const float Lambda = m_microsurfaceslope->Lambda(wi);
	// value
	const float value = 1.0f / (1.0f + Lambda);
	return value;	
}

float Microsurface::sampleHeight(const Vector3& wr, const float hr, const float U) const
{
	const float projectedArea = m_microsurfaceslope->projectedArea(-wr);

	if(projectedArea < 0.00001f)
		return (wr.z < 0.0f) ? hr : 0.0f;	

	const float dh = -logf(U) * std::abs(wr.z) / projectedArea;

	const float h = std::min(0.0f, hr) + dh * (wr.z>0.0f?1.0f:-1.0f);

	return h;
}

/* Microsurface made of conductor material */
class MicrosurfaceConductor : public Microsurface
{
public:
	const float m_n_i;
	const float m_n;
	const float m_k;
	MicrosurfaceConductor(
				 const bool slope_beckmann, // Beckmann or GGX
				 const float alpha_x,
				 const float alpha_y,
				 const float n_i, // ior of outside medium
				 const float n, // real part of conductor ior
				 const float k) // complex part of conductor ior
		: Microsurface(slope_beckmann, alpha_x, alpha_y),
		m_n_i(n_i), m_n(n), m_k(k)
	{}

public:
	// evaluate local phase function 
	//virtual float evalPhaseFunction(const Vector3& wi, const Vector3& wo) const;
	// sample local phase function
	virtual Vector3 samplePhaseFunction(const Vector3& wi) const; 
	virtual Vector3 samplePhaseFunction2(const Vector3& wi, float& sampleWeight ) const;
	Vector3 sample( const Vector3& wi, int& scatteringOrder, float& sampleWeight ) const;

	// evaluate BSDF limited to single scattering 
	// this is in average equivalent to eval(wi, wo, 1);
	//virtual float evalSingleScattering(const vec3& wi, const vec3& wo) const; 


protected:
	float Fresnel(const Vector3& wi, const Vector3& wm, const float n_i, const float n, const float k) const;
};

float MicrosurfaceConductor::Fresnel(const Vector3& wi, const Vector3& wm, const float ni, const float n, const float k) const
{
	// reflectance from a smooth conductor

	const float u = Dot(wi,wm);

	const float p = Sqrt(-Power(k,2) + Power(n,2) - Power(ni,2) + Power(ni,2)*Power(u,2) + 
     Sqrt(4*Power(k,2)*Power(n,2) + 
       Power(-Power(k,2) + Power(n,2) + Power(ni,2)*(-1 + Power(u,2)),2)))/Sqrt(2);

	const float q = Sqrt(Power(k,2) - Power(n,2) + Power(ni,2) - Power(ni,2)*Power(u,2) + 
     Sqrt(4*Power(k,2)*Power(n,2) + 
       Power(-Power(k,2) + Power(n,2) + Power(ni,2)*(-1 + Power(u,2)),2)))/Sqrt(2);

	const float rhoPerp = (Power(q,2) + Power(-p + ni*u,2))/(Power(q,2) + Power(p + ni*u,2));
	const float rhoPar = ((Power(q,2) + Power(-p + ni*u,2))*(Power(q,2) + Power(p - (ni*(1 - Power(u,2)))/u,2)))/
   ((Power(q,2) + Power(p + ni*u,2))*(Power(q,2) + Power(p + (ni*(1 - Power(u,2)))/u,2)));

	return 0.5f * ( rhoPerp + rhoPar );

}

Vector3 MicrosurfaceConductor::sample( const Vector3& wi, int& scatteringOrder, float& sampleWeight ) const
{
	// init
	Vector3 wr = -wi;
	float hr = 0.0f;
	
	// random walk
	scatteringOrder = 0;
	sampleWeight = 1.0f;
	while(true)
	{
		// next height
		float U = RandomReal();
		hr = sampleHeight(wr, hr, U);

		// leave the microsurface?
		if( hr >= 0.0f )
			break;
		else
			scatteringOrder++;

		// next direction
		wr = samplePhaseFunction2(-wr, sampleWeight);

		// if NaN (should not happen, just in case)
		if( (hr != hr) || (wr.z != wr.z) ) 
			return Vector3(0,0,1);
	}

	return wr;
}

Vector3 MicrosurfaceConductor::samplePhaseFunction(const Vector3& wi) const
{
	const float U1 = RandomReal();
	const float U2 = RandomReal();

	Vector3 wm = m_microsurfaceslope->sampleD_wi(wi, U1, U2);

	// reflect
	const Vector3 wo = -wi + 2.0f * wm * Dot(wi, wm);

	return wo;
}

Vector3 MicrosurfaceConductor::samplePhaseFunction2(const Vector3& wi, float& sampleWeight ) const
{
	const float U1 = RandomReal();
	const float U2 = RandomReal();

	Vector3 wm = m_microsurfaceslope->sampleD_wi(wi, U1, U2);
	sampleWeight *= Fresnel( wi, wm, m_n_i, m_n, m_k );

	// reflect
	const Vector3 wo = -wi + 2.0f * wm * Dot(wi, wm);

	return wo;
}

class MicrosurfaceMirror : public Microsurface
{
public:
	MicrosurfaceMirror(
				 const bool slope_beckmann, // Beckmann or GGX
				 const float alpha_x,
				 const float alpha_y ) 
		: Microsurface(slope_beckmann, alpha_x, alpha_y)
	{}

public:
	virtual Vector3 samplePhaseFunction(const Vector3& wi) const; 
	virtual Vector3 samplePhaseFunction2(const Vector3& wi, float& sampleWeight ) const;
	Vector3 sample( const Vector3& wi, int& scatteringOrder, float& sampleWeight ) const;

};

Vector3 MicrosurfaceMirror::samplePhaseFunction(const Vector3& wi) const
{
	const float U1 = RandomReal();
	const float U2 = RandomReal();

	Vector3 wm = m_microsurfaceslope->sampleD_wi(wi, U1, U2);

	// reflect
	const Vector3 wo = -wi + 2.0f * wm * Dot(wi, wm);

	return wo;
}

Vector3 MicrosurfaceMirror::samplePhaseFunction2(const Vector3& wi, float& sampleWeight ) const
{
	const float U1 = RandomReal();
	const float U2 = RandomReal();

	Vector3 wm = m_microsurfaceslope->sampleD_wi(wi, U1, U2);

	// reflect
	const Vector3 wo = -wi + 2.0f * wm * Dot(wi, wm);

	return wo;
}

Vector3 MicrosurfaceMirror::sample( const Vector3& wi, int& scatteringOrder, float& sampleWeight ) const
{
	// init
	Vector3 wr = -wi;
	float hr = 0.0f;
	
	// random walk
	scatteringOrder = 0;
	sampleWeight = 1.0f;
	while(true)
	{
		// next height
		float U = RandomReal();
		hr = sampleHeight(wr, hr, U);

		// leave the microsurface?
		if( hr >= 0.0f )
			break;
		else
			scatteringOrder++;

		// next direction
		wr = samplePhaseFunction2(-wr, sampleWeight);

		// if NaN (should not happen, just in case)
		if( (hr != hr) || (wr.z != wr.z) ) 
			return Vector3(0,0,1);
	}

	return wr;
}

class MicrosurfaceDielectric : public Microsurface
{
public:
	const float m_eta;
public:
	MicrosurfaceDielectric(const bool slope_beckmann, // Beckmann or GGX
				 const float alpha_x,
				 const float alpha_y,
				 const float eta)
		: Microsurface(slope_beckmann, alpha_x, alpha_y),
		m_eta(eta)
	{}

	// evaluate BSDF with a random walk (stochastic but unbiased)
	// scatteringOrder=0 --> contribution from all scattering events
	// scatteringOrder=1 --> contribution from 1st bounce only
	// scatteringOrder=2 --> contribution from 2nd bounce only, etc..
	//virtual float eval(const Vector3& wi, const Vector3& wo, const int scatteringOrder=0) const; 

	// evaluate BSDF with bidirectional algorithm for variance reduction
	//virtual float evalBidir(const Vector3& wi, const Vector3& wo, const int scatteringOrder=0) const;
	//virtual float evalMIS(const Vector3& wi, const Vector3& wo, const int scatteringOrder) const; 

	// sample final BSDF with a random walk
	// scatteringOrder is set to the number of bounces computed for this sample
	virtual Vector3 sample(const Vector3& wi, int& scatteringOrder) const;

public:
	// evaluate local phase function 
	virtual float evalPhaseFunction(const Vector3& wi, const Vector3& wo) const;
	float evalPhaseFunction(const Vector3& wi, const Vector3& wo, const bool wi_outside, const bool wo_outside) const;
	float MISweight(const Vector3& wi, const Vector3& wo, const bool wi_outside, const bool wo_outside) const;
	// sample local phase function
	virtual Vector3 samplePhaseFunction(const Vector3& wi) const; 
	Vector3 samplePhaseFunction(const Vector3& wi, const bool wi_outside, bool& wo_outside) const; 

	// evaluate BSDF limited to single scattering 
	// this is in average equivalent to eval(wi, wo, 1);
	//virtual float evalSingleScattering(const Vector3& wi, const Vector3& wo) const; 

protected:
	float Fresnel(const Vector3& wi, const Vector3& wm, const float eta) const;
	Vector3 refract(const Vector3 &wi, const Vector3 &wm, const float eta) const;
};

Vector3 MicrosurfaceDielectric::refract(const Vector3 &wi, const Vector3 &wm, const float eta) const 
{
	const float cos_theta_i = Dot(wi, wm);
	const float cos_theta_t2 = 1.0f - (1.0f-cos_theta_i*cos_theta_i) / (eta*eta);
	const float cos_theta_t = -sqrtf(std::max(0.0f,cos_theta_t2));

	return wm * (Dot(wi, wm) / eta + cos_theta_t) - wi / eta;
}

float MicrosurfaceDielectric::evalPhaseFunction(const Vector3& wi, const Vector3& wo) const
{
	return evalPhaseFunction(wi, wo, true, true) + evalPhaseFunction(wi, wo, true, false);
}

float MicrosurfaceDielectric::evalPhaseFunction(const Vector3& wi, const Vector3& wo, const bool wi_outside, const bool wo_outside) const
{
	const float eta = wi_outside ? m_eta : 1.0f / m_eta;

	if( wi_outside == wo_outside ) // reflection
	{
		// half vector 
		const Vector3 wh = Normalize(wi+wo);
		// value
		const float value = (wi_outside) ?
                       (0.25f * m_microsurfaceslope->D_wi(wi, wh) / Dot(wi, wh) * Fresnel(wi, wh, eta)) :
                       (0.25f * m_microsurfaceslope->D_wi(-wi, -wh) / Dot(-wi, -wh) * Fresnel(-wi, -wh, eta)) ;
		return value;
	}
	else // transmission
	{
		Vector3 wh = -Normalize(wi+wo*eta);
		wh = wh * ((wi_outside) ? double(sign(wh.z)) : double(-sign(wh.z)));

		if(Dot(wh, wi) < 0)
			return 0;

		float value;
		if(wi_outside){
			value = eta*eta * (1.0f-Fresnel(wi, wh, eta)) *
				m_microsurfaceslope->D_wi(wi, wh) * std::max(0.0, -Dot(wo, wh)) *
				1.0f / powf(Dot(wi, wh)+eta*Dot(wo,wh), 2.0f);
		}
		else
		{
			value = eta*eta * (1.0f-Fresnel(-wi, -wh, eta)) * 
				m_microsurfaceslope->D_wi(-wi, -wh) * std::max(0.0, -Dot(-wo, -wh)) *
				1.0f / powf(Dot(-wi, -wh)+eta*Dot(-wo,-wh), 2.0f);
		}

		return value;	
	}
}

float MicrosurfaceDielectric::MISweight(const Vector3& wi, const Vector3& wo, const bool wi_outside, const bool wo_outside) const
{
	const float eta = wi_outside ? m_eta : 1.0f / m_eta;

	if( wi_outside == wo_outside ) // reflection
	{
		// half vector 
		const Vector3 wh = Normalize(wi+wo);
		const float value =  m_microsurfaceslope->D(wh * sign(wh.z));
		return value;
	}
	else // transmission
	{
		Vector3 wh = Normalize(wi+wo*eta);
		const float value =  m_microsurfaceslope->D(wh * sign(wh.z));
		return value;	
	}
}


float MicrosurfaceDielectric::Fresnel(const Vector3& wi, const Vector3& wm, const float eta) const
{	
	const float cos_theta_i = Dot(wi, wm);
	const float cos_theta_t2 = 1.0f - (1.0f-cos_theta_i*cos_theta_i) / (eta*eta);

	// total internal reflection 
	if (cos_theta_t2 <= 0.0f) return 1.0f;

	const float cos_theta_t = sqrtf(cos_theta_t2);

	const float Rs = (cos_theta_i - eta * cos_theta_t) / (cos_theta_i + eta * cos_theta_t);
	const float Rp = (eta * cos_theta_i - cos_theta_t) / (eta * cos_theta_i + cos_theta_t);

	const float F = 0.5f * (Rs * Rs + Rp * Rp);
	return F;
}

Vector3 MicrosurfaceDielectric::samplePhaseFunction(const Vector3& wi) const
{
	bool wo_outside;
	return samplePhaseFunction(wi, true, wo_outside);
}

Vector3 MicrosurfaceDielectric::samplePhaseFunction(const Vector3& wi, const bool wi_outside, bool& wo_outside) const
{
	const float U1 = RandomReal();
	const float U2 = RandomReal();

	const float eta = wi_outside ? m_eta : 1.0f / m_eta;

	Vector3 wm = wi_outside ? (m_microsurfaceslope->sampleD_wi(wi, U1, U2)) :
						   (-m_microsurfaceslope->sampleD_wi(-wi, U1, U2)) ;

	const float F = Fresnel(wi, wm, eta);

	if( RandomReal() < F )
	{
		const Vector3 wo = -wi + 2.0f * wm * Dot(wi, wm); // reflect
		return wo;
	}
	else
	{
		wo_outside = !wi_outside;
		const Vector3 wo = refract(wi, wm, eta);
		return Normalize(wo);
	}
}

Vector3 MicrosurfaceDielectric::sample(const Vector3& wi, int& scatteringOrder) const
{
	// init
	Vector3 wr = -wi;
	float hr_outside = 0.0f;
	float hr_inside = -10000.0f; //////////////
	bool outside = true;
	
	// random walk
	scatteringOrder = 0;	
	while(true)
	{
		// next height
		float U = RandomReal();
		if(outside)
		{
			hr_outside = sampleHeight(wr, hr_outside, U);
			hr_inside = logf(1.0f - expf(hr_outside));
		}
		else
		{
			hr_inside = sampleHeight(-wr, hr_inside, U);
			hr_outside = logf(1.0f - expf(hr_inside));
		}

		// leave the microsurface?
		if( hr_outside >= 0.0f || hr_inside >= 0.0f)
			break;
		else
			scatteringOrder++;

		// next direction
		wr = samplePhaseFunction(-wr, outside, outside);		

		// if NaN (should not happen, just in case)
		if( (hr_outside != hr_outside) || (hr_inside != hr_inside)  || (wr.z != wr.z) ) 
			return Vector3(0,0,1);
	}

	return wr;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Diffuse
//////////////////////////////////////////////////////////////////////////////////////////////////
class MicrosurfaceDiffuse : public Microsurface
{
public:
	MicrosurfaceDiffuse(const bool slope_beckmann, // Beckmann or GGX
				 const float alpha_x,
				 const float alpha_y)
		: Microsurface(slope_beckmann, alpha_x, alpha_y)
	{}

public:
	// evaluate local phase function 
	//virtual float evalPhaseFunction(const Vector3& wi, const Vector3& wo) const;
	// sample local phase function
	virtual Vector3 samplePhaseFunction(const Vector3& wi) const; 

	// evaluate BSDF limited to single scattering 
	// this is in average equivalent to eval(wi, wo, 1);
	//virtual float evalSingleScattering(const Vector3& wi, const Vector3& wo) const; 
};


Vector3 MicrosurfaceDiffuse::samplePhaseFunction(const Vector3& wi) const
{
	const float U1 = RandomReal();
	const float U2 = RandomReal();
	const float U3 = RandomReal();
	const float U4 = RandomReal();

	Vector3 wm = m_microsurfaceslope->sampleD_wi(wi, U1, U2);

	// sample diffuse reflection
	Vector3 w1(0,0,0);
	Vector3 w2(0,0,0);
	buildOrthonormalBasis(w1, w2, wm);

	float r1 = 2.0f*U3 - 1.0f;
	float r2 = 2.0f*U4 - 1.0f;

	// concentric map code from
	// http://psgraphics.blogspot.ch/2011/01/improved-code-for-concentric-map.html
	float phi, r;
	if (r1 == 0 && r2 == 0) {
		r = phi = 0;
	} else if (r1*r1 > r2*r2) {
		r = r1;
		phi = (M_PI/4.0f) * (r2/r1);
	} else {
		r = r2;
		phi = (M_PI/2.0f) - (r1/r2) * (M_PI/4.0f);
	}
	float x = r*cosf(phi);
	float y = r*sinf(phi);
	float z = sqrtf(std::max(0.0f, 1.0f - x*x - y*y));
	Vector3 wo = x*w1 + y*w2 + z*wm;

	return wo;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Biscale - A microsurface where each microfacet has a rough mirror BRDF
//  microfacets do not support anisotropy
//////////////////////////////////////////////////////////////////////////////////////////////////
class MicrosurfaceBiScaleMirror : public Microsurface
{
public:
	MicrosurfaceBiScaleMirror(
														const bool macro_slope_beckmann, // Beckmann or GGX
														const float macro_alpha_x,
														const float macro_alpha_y,
														const bool micro_slope_beckmann, // Beckmann or GGX
														const float micro_alpha
														)
		: Microsurface(macro_slope_beckmann, macro_alpha_x, macro_alpha_y), 
		  m_microfacetBRDF( micro_slope_beckmann, micro_alpha, micro_alpha )
	{}

public:
	virtual Vector3 samplePhaseFunction(const Vector3& wi) const;

private:
	MicrosurfaceMirror m_microfacetBRDF;
};

Vector3 MicrosurfaceBiScaleMirror::samplePhaseFunction(const Vector3& wi) const
{
	const float U1 = RandomReal();
	const float U2 = RandomReal();

	const Vector3 wm = m_microsurfaceslope->sampleD_wi(wi, U1, U2);

	// sample BRDF
	Vector3 w1(0,0,0);
	Vector3 w2(0,0,0);
	buildOrthonormalBasis(w1, w2, wm);

	int order(0);
	float weight(0.0);
	// rotate wi into frame of microfacet and call microfacet.sample()
	const Vector3 dir = m_microfacetBRDF.sample( Vector3( Dot( wi, w1 ), Dot( wi, w2 ), Dot( wi, wm ) ), order, weight );

	// rotate the sampled direction back
	const Vector3 wo = dir.x * w1 + dir.y * w2 + dir.z * wm;
	return wo;
}