using UnityEngine;
using System.Collections;
using System.Collections.Generic;

public class fts
{
    private static bool IsZero(double d)
    {
        const double eps = 1e-9;
        return d > -eps && d < eps;
    }

    private static double GetCubicRoot(double value)
    {
        if (value > 0.0)
        {
            return System.Math.Pow(value, 1.0 / 3.0);
        }
        else if (value < 0)
        {
            return -System.Math.Pow(-value, 1.0 / 3.0);
        }
        else
        {
            return 0.0;
        }
    }

    public static int SolveQuadric(double c0, double c1, double c2, out double s0, out double s1)
    {
        s0 = double.NaN;
        s1 = double.NaN;

        double p, q, D;

        p = c1 / (2 * c0);
        q = c2 / c0;

        D = p * p - q;

        if (IsZero(D))
        {
            s0 = -p;
            return 1;
        }
        else if (D < 0)
        {
            return 0;
        }
        else
        {
            double sqrt_D = System.Math.Sqrt(D);

            s0 = sqrt_D - p;
            s1 = -sqrt_D - p;
            return 2;
        }
    }

    public static int SolveCubic(double c0, double c1, double c2, double c3, out double s0, out double s1, out double s2)
    {
        s0 = double.NaN;
        s1 = double.NaN;
        s2 = double.NaN;

        int num;
        double sub;
        double A, B, C;
        double sq_A, p, q;
        double cb_p, D;

        A = c1 / c0;
        B = c2 / c0;
        C = c3 / c0;

        sq_A = A * A;
        p = 1.0 / 3 * (-1.0 / 3 * sq_A + B);
        q = 1.0 / 2 * (2.0 / 27 * A * sq_A - 1.0 / 3 * A * B + C);

        cb_p = p * p * p;
        D = q * q + cb_p;

        if (IsZero(D))
        {
            if (IsZero(q))
            {
                s0 = 0;
                num = 1;
            }
            else
            {
                double u = GetCubicRoot(-q);
                s0 = 2 * u;
                s1 = -u;
                num = 2;
            }
        }
        else if (D < 0)
        {
            double phi = 1.0 / 3 * System.Math.Acos(-q / System.Math.Sqrt(-cb_p));
            double t = 2 * System.Math.Sqrt(-p);

            s0 = t * System.Math.Cos(phi);
            s1 = -t * System.Math.Cos(phi + System.Math.PI / 3);
            s2 = -t * System.Math.Cos(phi - System.Math.PI / 3);
            num = 3;
        }
        else
        {
            double sqrt_D = System.Math.Sqrt(D);
            double u = GetCubicRoot(sqrt_D - q);
            double v = -GetCubicRoot(sqrt_D + q);

            s0 = u + v;
            num = 1;
        }

        sub = 1.0 / 3 * A;

        if (num > 0) s0 -= sub;
        if (num > 1) s1 -= sub;
        if (num > 2) s2 -= sub;

        return num;
    }

    public static int SolveQuartic(double c0, double c1, double c2, double c3, double c4, out double s0, out double s1, out double s2, out double s3)
    {
        s0 = double.NaN;
        s1 = double.NaN;
        s2 = double.NaN;
        s3 = double.NaN;

        double[] coeffs = new double[4];
        double z, u, v, sub;
        double A, B, C, D;
        double sq_A, p, q, r;
        int num;

        A = c1 / c0;
        B = c2 / c0;
        C = c3 / c0;
        D = c4 / c0;

        sq_A = A * A;
        p = -3.0 / 8 * sq_A + B;
        q = 1.0 / 8 * sq_A * A - 1.0 / 2 * A * B + C;
        r = -3.0 / 256 * sq_A * sq_A + 1.0 / 16 * sq_A * B - 1.0 / 4 * A * C + D;

        if (IsZero(r))
        {
            coeffs[3] = q;
            coeffs[2] = p;
            coeffs[1] = 0;
            coeffs[0] = 1;

            num = fts.SolveCubic(coeffs[0], coeffs[1], coeffs[2], coeffs[3], out s0, out s1, out s2);
        }
        else
        {
            coeffs[3] = 1.0 / 2 * r * p - 1.0 / 8 * q * q;
            coeffs[2] = -r;
            coeffs[1] = -1.0 / 2 * p;
            coeffs[0] = 1;

            SolveCubic(coeffs[0], coeffs[1], coeffs[2], coeffs[3], out s0, out s1, out s2);

            z = s0;

            u = z * z - r;
            v = 2 * z - p;

            if (IsZero(u))
                u = 0;
            else if (u > 0)
                u = System.Math.Sqrt(u);
            else
                return 0;

            if (IsZero(v))
                v = 0;
            else if (v > 0)
                v = System.Math.Sqrt(v);
            else
                return 0;

            coeffs[2] = z - u;
            coeffs[1] = q < 0 ? -v : v;
            coeffs[0] = 1;

            num = fts.SolveQuadric(coeffs[0], coeffs[1], coeffs[2], out s0, out s1);

            coeffs[2] = z + u;
            coeffs[1] = q < 0 ? v : -v;
            coeffs[0] = 1;

            if (num == 0) num += fts.SolveQuadric(coeffs[0], coeffs[1], coeffs[2], out s0, out s1);
            else if (num == 1) num += fts.SolveQuadric(coeffs[0], coeffs[1], coeffs[2], out s1, out s2);
            else if (num == 2) num += fts.SolveQuadric(coeffs[0], coeffs[1], coeffs[2], out s2, out s3);
        }

        sub = 1.0 / 4 * A;

        if (num > 0) s0 -= sub;
        if (num > 1) s1 -= sub;
        if (num > 2) s2 -= sub;
        if (num > 3) s3 -= sub;

        return num;
    }

    public static float ballistic_range(float speed, float gravity, float initial_height)
    {
        Debug.Assert(speed > 0 && gravity > 0 && initial_height >= 0, "fts.ballistic_range called with invalid data");

        float angle = 45 * Mathf.Deg2Rad;
        float cos = Mathf.Cos(angle);
        float sin = Mathf.Sin(angle);

        float range = (speed * cos / gravity) * (speed * sin + Mathf.Sqrt(speed * speed * sin * sin + 2 * gravity * initial_height));
        return range;
    }

    public static int solve_ballistic_arc(Vector3 proj_pos, float proj_speed, Vector3 target, float gravity, out Vector3 s0, out Vector3 s1)
    {
        Debug.Assert(proj_pos != target && proj_speed > 0 && gravity > 0, "fts.solve_ballistic_arc called with invalid data");

        s0 = Vector3.zero;
        s1 = Vector3.zero;

        Vector3 diff = target - proj_pos;
        Vector3 diffXZ = new Vector3(diff.x, 0f, diff.z);
        float groundDist = diffXZ.magnitude;

        float speed2 = proj_speed * proj_speed;
        float speed4 = proj_speed * proj_speed * proj_speed * proj_speed;
        float y = diff.y;
        float x = groundDist;
        float gx = gravity * x;

        float root = speed4 - gravity * (gravity * x * x + 2 * y * speed2);

        if (root < 0)
            return 0;

        root = Mathf.Sqrt(root);

        float lowAng = Mathf.Atan2(speed2 - root, gx);
        float highAng = Mathf.Atan2(speed2 + root, gx);
        int numSolutions = lowAng != highAng ? 2 : 1;

        Vector3 groundDir = diffXZ.normalized;
        s0 = groundDir * Mathf.Cos(lowAng) * proj_speed + Vector3.up * Mathf.Sin(lowAng) * proj_speed;
        if (numSolutions > 1)
            s1 = groundDir * Mathf.Cos(highAng) * proj_speed + Vector3.up * Mathf.Sin(highAng) * proj_speed;

        return numSolutions;
    }

    public static int solve_ballistic_arc(Vector3 proj_pos, float proj_speed, Vector3 target_pos, Vector3 target_velocity, float gravity, out Vector3 s0, out Vector3 s1)
    {
        s0 = Vector3.zero;
        s1 = Vector3.zero;

        double G = gravity;

        double A = proj_pos.x;
        double B = proj_pos.y;
        double C = proj_pos.z;
        double M = target_pos.x;
        double N = target_pos.y;
        double O = target_pos.z;
        double P = target_velocity.x;
        double Q = target_velocity.y;
        double R = target_velocity.z;
        double S = proj_speed;

        double H = M - A;
        double J = O - C;
        double K = N - B;
        double L = -.5f * G;

        double c0 = L * L;
        double c1 = -2 * Q * L;
        double c2 = Q * Q - 2 * K * L - S * S + P * P + R * R;
        double c3 = 2 * K * Q + 2 * H * P + 2 * J * R;
        double c4 = K * K + H * H + J * J;

        double[] times = new double[4];
        int numTimes = SolveQuartic(c0, c1, c2, c3, c4, out times[0], out times[1], out times[2], out times[3]);

        System.Array.Sort(times);

        Vector3[] solutions = new Vector3[2];
        int numSolutions = 0;

        for (int i = 0; i < times.Length && numSolutions < 2; ++i)
        {
            double t = times[i];
            if (t <= 0 || double.IsNaN(t))
                continue;

            solutions[numSolutions].x = (float)((H + P * t) / t);
            solutions[numSolutions].y = (float)((K + Q * t - L * t * t) / t);
            solutions[numSolutions].z = (float)((J + R * t) / t);
            ++numSolutions;
        }

        if (numSolutions > 0) s0 = solutions[0];
        if (numSolutions > 1) s1 = solutions[1];

        return numSolutions;
    }

    public static bool solve_ballistic_arc_lateral(Vector3 proj_pos, float lateral_speed, Vector3 target_pos, float max_height, out Vector3 fire_velocity, out float gravity)
    {
        Debug.Assert(proj_pos != target_pos && lateral_speed > 0 && max_height > proj_pos.y, "fts.solve_ballistic_arc called with invalid data");

        fire_velocity = Vector3.zero;
        gravity = float.NaN;

        Vector3 diff = target_pos - proj_pos;
        Vector3 diffXZ = new Vector3(diff.x, 0f, diff.z);
        float lateralDist = diffXZ.magnitude;

        if (lateralDist == 0)
            return false;

        float time = lateralDist / lateral_speed;

        fire_velocity = diffXZ.normalized * lateral_speed;

        float a = proj_pos.y;
        float b = max_height;
        float c = target_pos.y;

        gravity = -4 * (a - 2 * b + c) / (time * time);
        fire_velocity.y = -(3 * a - 4 * b + c) / time;

        return true;
    }

    public static bool solve_ballistic_arc_lateral(Vector3 proj_pos, float lateral_speed, Vector3 target, Vector3 target_velocity, float max_height_offset, out Vector3 fire_velocity, out float gravity, out Vector3 impact_point)
    {
        Debug.Assert(proj_pos != target && lateral_speed > 0, "fts.solve_ballistic_arc_lateral called with invalid data");

        fire_velocity = Vector3.zero;
        gravity = 0f;
        impact_point = Vector3.zero;

        Vector3 targetVelXZ = new Vector3(target_velocity.x, 0f, target_velocity.z);
        Vector3 diffXZ = target - proj_pos;
        diffXZ.y = 0;

        float c0 = Vector3.Dot(targetVelXZ, targetVelXZ) - lateral_speed * lateral_speed;
        float c1 = 2f * Vector3.Dot(diffXZ, targetVelXZ);
        float c2 = Vector3.Dot(diffXZ, diffXZ);
        double t0, t1;
        int n = fts.SolveQuadric(c0, c1, c2, out t0, out t1);

        bool valid0 = n > 0 && t0 > 0;
        bool valid1 = n > 1 && t1 > 0;

        float t;
        if (!valid0 && !valid1)
            return false;
        else if (valid0 && valid1)
            t = Mathf.Min((float)t0, (float)t1);
        else
            t = valid0 ? (float)t0 : (float)t1;

        impact_point = target + (target_velocity * t);

        Vector3 dir = impact_point - proj_pos;
        fire_velocity = new Vector3(dir.x, 0f, dir.z).normalized * lateral_speed;

        float a = proj_pos.y; 
        float b = Mathf.Max(proj_pos.y, impact_point.y) + max_height_offset;
        float c = impact_point.y;

        gravity = -4 * (a - 2 * b + c) / (t * t);
        fire_velocity.y = -(3 * a - 4 * b + c) / t;

        return true;
    }
}