using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using static FunctionLibrary;
using static UnityEngine.Mathf;

public static class FunctionLibrary
{
    public delegate Vector3 Function(float u, float v, float t);

    public enum FunctionName { MultiWave, Wave, Ripple, Sphere, MorphedSphere, EntangledSphere, Torus, Diffusion, Helix, Helicoid, EllipticalHelicoid, Lemniscate, FractalArt, Paraboloid };

    public static Function[] functions = { MultiWave, Wave, Ripple, Sphere, MorphedSphere, EntangledSphere, Torus, Diffusion, Helix, Helicoid, EllipticalHelicoid, Lemniscate, FractalArt, Paraboloid };

    public static FunctionName GetNextFunctionName(FunctionName name)
    {
        return (int)name < functions.Length - 1 ? name + 1 : 0;
    }

    public static Function GetFunction(FunctionName name)
    {
        return functions[(int)name];
    }

    public static FunctionName GetRandomFunctionNameOtherThan(FunctionName name)
    {
        var choice = (FunctionName)Random.Range(1, functions.Length);
        return choice == name ? 0 : choice;
    }

    public static Vector3 Morph(float u, float v, float t, Function from, Function to, float progress)
    {
        return Vector3.LerpUnclamped(from(u, v, t), to(u, v, t), SmoothStep(0f, 1f, progress));
    }

    public static Vector3 MultiWave(float u, float v, float t)
    {
        Vector3 p;
        p.x = u;
        p.y = Sin(PI * (u + 0.5f * t));
        p.y += 0.5f * Sin(2f * PI * (v + t));
        p.y += Sin(PI * (u + v + 0.25f * t));
        p.y *= 1f / 2.5f;
        p.z = v;
        return p;
    }

    public static Vector3 Wave(float u, float v, float t)
    {
        Vector3 p;
        p.x = u;
        p.y = Sin(PI * (u + v + t));
        p.z = v;
        return p;
    }

    public static Vector3 Ripple(float u, float v, float t)
    {
        float d = Sqrt(u * u + v * v);
        Vector3 p;
        p.x = u;
        p.y = Sin(PI * (4f * d - t));
        p.y /= 1f + 10f * d;
        p.z = v;
        return p;
    }

    public static Vector3 Sphere(float u, float v, float t)
    {
        Vector3 p;
        float radius = 1f;
        float amplitude = 0.3f;
        float speed = 2f;
        float theta = PI * 2f * u;
        float phi = PI * (v - 0.5f);
        float r = radius + amplitude * Sin(speed * t);
        p.x = r * Cos(theta) * Cos(phi);
        p.y = r * Sin(phi);
        p.z = r * Sin(theta) * Cos(phi);
        return p;
    }

    public static Vector3 MorphedSphere(float u, float v, float t)
    {
        float r = 0.9f + 0.1f * Sin(PI * (6f * u + 4f * v + t));
        float s = r * Cos(0.5f * PI * v);
        Vector3 p;
        p.x = s * Sin(PI * u);
        p.y = r * Sin(PI * 0.5f * v);
        p.z = s * Cos(PI * u);
        return p;
    }

    public static Vector3 EntangledSphere(float u, float v, float t)
    {
        Vector3 p;
        float a = 2f * PI * u;
        float b = 2f * PI * v;
        float s = 0.2f * Sin(PI * t);
        float r = 1 - s * Cos(a) * Sin(b) - s / 2f * Sin(a);
        p.x = r * Sin(a) * Sin(b);
        p.y = r * Cos(b);
        p.z = r * Cos(a) * Sin(b);
        return p;
    }

    public static Vector3 Torus(float u, float v, float t)
    {
        float r1 = 0.7f + 0.1f * Sin(PI * (6f * u + 0.5f *t));
        float r2 = 0.15f + 0.05f * Sin(PI * 8f * u + 4f * v +2f * t);
        float s = r1 + r2 * Cos(PI * v);
        Vector3 p;
        p.x = s * Sin(PI * u);
        p.y = r2 * Sin(PI * v);
        p.z = s * Cos(PI * u);
        return p;
    }

    public static Vector3 Diffusion(float u, float v, float t)
    {
        Vector3 p;
        float s = 0.1f + Exp(-4 * (u * u + v * v) / (0.001f + t % 1.0f));
        p.x = u;
        p.y = s;
        p.z = v;
        return p;
    }

    public static Vector3 Helix(float u, float v, float t)
    {
        Vector3 p;
        float r = 0.3f + 0.1f * Sin(PI * (6f * v + t));
        float a = PI * (u - 0.5f * t);
        p.x = r * Sin(a);
        p.y = v;
        p.z = r * Cos(a);
        return p;
    }

    public static Vector3 Helicoid(float u, float v, float t)
    {
        Vector3 p;
        float alpha = PI * t;
        float cosAlpha = Cos(alpha);
        float sinAlpha = Sin(alpha);
        float sinhV = 0.5f * (Exp(v) - Exp(-v));
        float coshV = 0.5f * (Exp(v) + Exp(-v));
        float sinU = Sin(u);
        float cosU = Cos(u);
        p.x = cosAlpha * sinhV * sinU + sinAlpha * coshV * cosU;
        p.y = -cosAlpha * sinhV * cosU + sinAlpha * coshV * sinU;
        p.z = u * cosAlpha + v * sinAlpha;
        return p;
    }

    public static Vector3 EllipticalHelicoid(float u, float v, float t)
    {
        Vector3 p;
        float ut = u + t;
        p.x = v * Cos(PI * ut);
        p.y = v * Sin(PI * ut);
        p.z = u;
        return p;
    }

    public static Vector3 Lemniscate(float u, float v, float t)
    {
        Vector3 p;
        float a = PI * (u - t);
        float denom = Sqrt(2) * Pow(Sin(a), 2) + 1;
        p.x = Cos(a) / denom;
        p.y = Sin(a) * Cos(a) / denom;
        p.z = v;
        return p;
    }

    public static Vector3 FractalArt(float u, float v, float t)
    {
        Vector3 p;
        float scale = 0.08f;
        float frequency = 2f;
        int iterations = 8;
        Vector3 Fractal(Vector3 point, int depth)
        {
            if (depth <= 0)
                return point;
            point.x += scale * Cos(frequency * point.y + t);
            point.y += scale * Sin(frequency * point.x + t);
            point.z += scale * Sin(frequency * point.x * point.y + t);
            return Fractal(point, depth - 1);
        }
        Vector3 startPoint = new Vector3(u, v, 0f);
        p = Fractal(startPoint, iterations);
        return p;
    }

    public static Vector3 Paraboloid(float u, float v, float t)
    {
        Vector3 p;
        float scale = 1f;
        float height = 1f;
        float waveAmplitude = 0.2f;
        float waveFrequency = 3f;
        float x = u;
        float y = v * height;
        float z = (u * u - v * v) + waveAmplitude * Sin(waveFrequency);
        p.x = x * scale;
        p.y = y * scale;
        p.z = z * scale;
        return p;
    }
}
