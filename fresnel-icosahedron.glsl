// shadertoy : https://www.shadertoy.com/view/XcSBRW

/*
Copyright 2024 Taylor Wampler

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), 
to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, 
and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

// things to edit
#define CAMERA_DISTANCE 5.0
#define MAX_COLLISIONS 5
#define SCENE_UP vec3(0.0, 0.0, 1.0)
#define ICO_ABSORPTION vec3(2.5, 1.0, 1.0)
#define ICO_SURFACE_FIRST_REFLECTION_COLOR vec3(1.0) 
#define ICO_INDEX 2.1
#define AIR_INDEX 1.0003
#define SQRT_SSA_SAMPLES 2

// Give a little nudge to reflected and refracted rays.
// I had to play around with this until the artifacts went away.
#define RAY_NUDGE 10e-5

// Tolerance for the triangle intersection test.
#define INTERSECTION_TOL 10e-5

// geometric constants
#define ICO_HEIGHT (sqrt(5.0) / 5.0)
#define PENT_R (2.0 * sqrt(5.0) / 5.0)

struct Camera       { vec3 eye, target; };
struct Ray          { vec3 point, direction; };
struct Triangle     { vec3 p1, p2, p3; };
struct Intersection { float t; vec3 normal; };

// Think of the icosahedron as four shifted xz plane slices. 
// A and L are in the first and last slices respectively.
// BCDEF and GHIJK are the pentagons, 
// they are reflections of each other about the x-axis in the plane of the slice.
// Also note that the first and second slices are reflections 
// of the last and third slices about the xz plane.
//
// 
// front view       back view
//       B            J       I
//  F         C           L    
//       A           K         H
//   E       D            G        
//
// https://mathworld.wolfram.com/RegularIcosahedron.html
const vec3 pA = vec3(0.0,
               -1.0, 
               0.0);
const vec3 pB = vec3(0.0, 
                -ICO_HEIGHT, 
                PENT_R);
const vec3 pC = vec3(PENT_R*sqrt(0.125*(5.0+sqrt(5.0))), 
                -ICO_HEIGHT, 
                PENT_R*0.25*(sqrt(5.0)-1.0));
const vec3 pD = vec3(PENT_R*sqrt(0.125*(5.0-sqrt(5.0))), 
                -ICO_HEIGHT, 
                PENT_R*0.25*(-sqrt(5.0)-1.0));
const vec3 pE = vec3(-PENT_R*sqrt(0.125*(5.0-sqrt(5.0))), 
                -ICO_HEIGHT, 
                PENT_R*0.25*(-sqrt(5.0)-1.0));
const vec3 pF = vec3(-PENT_R*sqrt(0.125*(5.0+sqrt(5.0))), 
                -ICO_HEIGHT, 
                PENT_R*0.25*(sqrt(5.0)-1.0));
const vec3 pG = vec3(0.0, 
                ICO_HEIGHT, 
                -PENT_R);
const vec3 pH = vec3(PENT_R*sqrt(0.125*(5.0+sqrt(5.0))), 
                ICO_HEIGHT, 
                PENT_R*0.25*(1.0-sqrt(5.0)));
const vec3 pI = vec3(PENT_R*sqrt(0.125*(5.0-sqrt(5.0))), 
                ICO_HEIGHT, 
                PENT_R*0.25*(sqrt(5.0)+1.0));
const vec3 pJ = vec3(-PENT_R*sqrt(0.125*(5.0-sqrt(5.0))), 
                ICO_HEIGHT, 
                PENT_R*0.25*(sqrt(5.0)+1.0));
const vec3 pK = vec3(-PENT_R*sqrt(0.125*(5.0+sqrt(5.0))), 
                ICO_HEIGHT, 
                PENT_R*0.25*(1.0-sqrt(5.0)));
const vec3 pL = vec3(0.0,
               1.0, 
               0.0);
Triangle triangles[20];

//---------------------------------------------------------------

Ray cameraRay(in Camera c, in vec2 uv)
{
    vec3 forward = normalize(c.target - c.eye);
    vec3 left    = normalize(cross(forward, SCENE_UP));
    vec3 up      = cross(left, forward);
    
    vec3 dy = (iResolution.y / iResolution.x) * (uv.y - 0.5) * up ;
    vec3 dx = (uv.x - 0.5) * left;
    return Ray(c.eye, normalize(dy - dx + forward));
}

Intersection closerIntersection(in Intersection i1, in Intersection i2)
{
   if (i2.t >= 0.0 && (i1.t < 0.0 || i2.t < i1.t)) return i2;
   else return i1;
}

// Real‐Time Rendering, 4th Edition by Tomas Akenine-Möller, Eric Haines, Naty Hoffman, 
// Angelo Pesce, Michał Iwanicki, and Sébastien Hillaire
Intersection rayIntersectTriangle(in Ray ry, in Triangle tr)
{
    vec3 e1 = tr.p3 - tr.p1;
    vec3 e2 = tr.p2 - tr.p1;
    vec3 s = ry.point - tr.p1;
    vec3 m = cross(s, ry.direction);
    vec3 normal = cross(e1, e2);
    
    // normal
    float dp = -dot(ry.direction, normal);
    // check if ray is parallel to the plane
    if (dp > -INTERSECTION_TOL && dp < INTERSECTION_TOL) 
        return Intersection(-1.0, vec3(0.0));
    
    float inv_dp = 1.0 / dp;
    vec3 coord = inv_dp * vec3(dot(normal, s), dot(m, e2), -dot(m, e1));
    if (coord.x < 0.0 || coord.y < 0.0 || coord.z < 0.0 || (coord.y + coord.z) > 1.0)
        return Intersection(-1.0, vec3(0.0));
        
    // The normal returned should be a unit vector anti-aligned with the ray.
    if (dp < 0.0) normal = -normal;
    normal = normalize(normal);
    
    return Intersection(coord.x, normal);    
}

Intersection rayIntersectIcosahedron(in Ray ry)
{
    Intersection intersection;
    intersection.t = -1.0;
    intersection.normal = vec3(0.0);
    for (int i = 0; i < 20; i++)
    {
        Intersection collide = rayIntersectTriangle(ry, triangles[i]);
        intersection = closerIntersection(intersection, collide);
    }
    
    return intersection;
}

// Use Schlick's model to approximate Fresnel reflectance.
//
// Real‐Time Rendering, 4th Edition by Tomas Akenine-Möller, Eric Haines, Naty Hoffman, 
// Angelo Pesce, Michał Iwanicki, and Sébastien Hillaire
float reflectance(in float n1, in float n2, in Ray ry, in Intersection intersection)
{
    float cosTheta = clamp(-dot(ry.direction, intersection.normal), 0.0, 1.0);
    
    // check for total internal reflection
    float sineThetaSquared = 1.0 - cosTheta * cosTheta;
    float n1n2Squared = n1 / n2;
    n1n2Squared *= n1n2Squared;
    if (sineThetaSquared * n1n2Squared > 1.0) return 1.0;
    
    float r0 = (n1 - n2) / (n1 + n2);
    r0 *= r0;
    float f = 1.0 - cosTheta;
    
    return r0 + (1.0 - r0) * f * f * f * f * f;
}

vec3 getEnvironmentColor(in Ray ry)
{
    return texture(iChannel0, vec3(ry.direction.x, ry.direction.z, -ry.direction.y)).rgb;
}

vec3 getColor(Ray ry)
{
    vec3 color = vec3(0.0);
    float multiplier = 1.0;
    float totInternalPath = 0.0;
    for (int i = 0; i < MAX_COLLISIONS; i++)
    {
	Intersection intersection = rayIntersectIcosahedron(ry);
        if (intersection.t < 0.0 || intersection.t >= 1e2)
        {
            if (i == 0) return getEnvironmentColor(ry);
            break;
        }
        
        // Move along the ray to the intersection point.
        ry.point += ry.direction * intersection.t;
        
        // Determine the appropriate interface. The skybox and icosahedron are the only objects.
        float n1, n2;
        if (i > 0)
        {
            n1 = ICO_INDEX;
            n2 = AIR_INDEX;
        }
        else
        {
            n1 = AIR_INDEX;
            n2 = ICO_INDEX;
        }
        
        // Fresnel reflectance and the Law of Reflection
        float r = reflectance(n1, n2, ry, intersection);
        vec3 reflectDirection = normalize(reflect(ry.direction, intersection.normal));
        Ray ryReflect = Ray(ry.point + reflectDirection * RAY_NUDGE, reflectDirection);
        // Fresnel transmittance and the Law of Refraction (Snell's Law)
        float m = 1.0 - r;
        vec3 refractDirection = normalize(refract(ry.direction, intersection.normal, n1 / n2));
        Ray ryRefract = Ray(ry.point + refractDirection * RAY_NUDGE, refractDirection);
        
        // On the first intersection, add the reflection contribution and follow the refracted ray.
        // On later intersections, add the refraction contribution and follow the internally reflected ray.
        if (i > 0)
        {
            // Bouguer–Lambert / Beer-Lambert Law
            totInternalPath += intersection.t;
            vec3 absorption = exp(-ICO_ABSORPTION * totInternalPath);
            
            // refraction contribution from outside the icosahedron
            vec3 environmentColor = getEnvironmentColor(ryRefract);
            color += (m * multiplier * absorption) * environmentColor;
            
            ry = ryReflect;
            multiplier *= r;
        }
        else
        {
            // reflection contribution from outside the icosahedron
            vec3 environmentColor = getEnvironmentColor(ryReflect);
            color += (r * ICO_SURFACE_FIRST_REFLECTION_COLOR) * environmentColor;
            
            ry = ryRefract;
            multiplier *= m;
        }
    }
    
    return color;
}

//---------------------------------------------------------------

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = fragCoord/iResolution.xy;
    
    // third-person camera
    float azimuthal = 3.14 / 2.0 + iTime / 3.14;
    float polar     = 3.14 / 2.0;
    Camera camera;
    camera.eye    = CAMERA_DISTANCE * vec3(sin(polar) * cos(azimuthal), sin(polar) * sin(azimuthal), cos(polar));
    camera.target = vec3(0.0, 0.0, 0.0);

    // triangles connecting the first and second slices
    triangles[0]  = Triangle(pA, pB, pC);
    triangles[1]  = Triangle(pA, pC, pD);
    triangles[2]  = Triangle(pA, pD, pE);
    triangles[3]  = Triangle(pA, pE, pF);
    triangles[4]  = Triangle(pA, pF, pB);
    // triangles connecting the second and third slices
    triangles[5]  = Triangle(pB, pJ, pI);
    triangles[6]  = Triangle(pB, pI, pC);
    triangles[7]  = Triangle(pC, pI, pH);
    triangles[8]  = Triangle(pC, pH, pD);
    triangles[9]  = Triangle(pD, pH, pG);
    triangles[10] = Triangle(pD, pG, pE);
    triangles[11] = Triangle(pE, pG, pK);
    triangles[12] = Triangle(pE, pK, pF);
    triangles[13] = Triangle(pF, pK, pJ);
    triangles[14] = Triangle(pF, pJ, pB);
    // triangles connecting the third and last slices
    triangles[15] = Triangle(pL, pJ, pI);
    triangles[16] = Triangle(pL, pI, pH);
    triangles[17] = Triangle(pL, pH, pG);
    triangles[18] = Triangle(pL, pG, pK);
    triangles[19] = Triangle(pL, pK, pJ);

    // A basic SSA scheme. Shadertoy doesn't automatically antialias. 
    vec3 color = vec3(0.0);
    for (int i = 0; i < SQRT_SSA_SAMPLES; i++)
    {
        for (int j = 0; j < SQRT_SSA_SAMPLES; j++) 
        {
            // Calculate the color from a subpixel ray.
            vec2 offset = vec2(i, j) / (iResolution.xy * float(SQRT_SSA_SAMPLES));
            color += getColor(cameraRay(camera, uv + offset));
         }
    }
    color /= float(SQRT_SSA_SAMPLES * SQRT_SSA_SAMPLES);
    
    // gamma correction
    color = pow(max(color, 0.0), vec3(1.0/2.2));
    
    fragColor = vec4(color, 1.0);
}
