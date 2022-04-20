using System;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Factorization;
using PolyLib;


namespace AstrolibNS
{
    public static class Astrolib
    {
        // Constants
        public static double pi = 4 * Math.Atan(1);             
        public static double a = 6378.137;                // [km] earth equatorial radius
        public static double e = 8.1819190842622e-2;	  // [] 
        public static double f = 1 / 298.257223563;       // [] earth flattening factor
        public static double mu = 3.986004418e5;          // [km^3/s^2] earth standard gravitational parameter

        // Utility functions
        public static double Rad2Deg(double theta)          // radians to degrees conversion
        {
            return theta * (180 / pi);
        }

        public static double Deg2Rad(double theta)          // degrees to radians conversion
        {
            return theta * (pi / 180);
        }

        public static Vector<double> GetGroundSiteMat(Vector<double> llaMat, double jd)     // cartisian coordinates of ground site from geodetic coordinates
        {   
            var V = Vector<double>.Build;
            var rSite = V.Dense(3);

            double lat = llaMat[0];
            double lon = llaMat[1];
            double alt = llaMat[2];
        
            double LSTrad = Astrolib.LocalSiderialTime(jd, lon);            // local siderial time (angle) from j2k time
        
            rSite[0] = ((a / Math.Sqrt(1 - (2 * f - Math.Pow(f, 2)) * Math.Pow(Math.Sin(lat), 2))) + alt) * Math.Cos(lat) * Math.Cos(LSTrad);
            rSite[1] = ((a / Math.Sqrt(1 - (2 * f - Math.Pow(f, 2)) * Math.Pow(Math.Sin(lat), 2))) + alt) * Math.Cos(lat) * Math.Sin(LSTrad);
            rSite[2] = (((a * Math.Pow(1-f, 2)) / Math.Sqrt(1 - (2 * f - Math.Pow(f, 2)) * Math.Pow(Math.Sin(lat), 2))) + alt) * Math.Sin(LSTrad);
            
            return rSite;
        }

        public static Vector<double> GetLOSMat(Vector<double> RA_Dec)  // Compute line of sight matrix from right ascension and declination
        {
            var V = Vector<double>.Build; 
            var LOS = V.Dense(3);

            double alpha = RA_Dec[0];
            double delta = RA_Dec[1];
            
            LOS[0] = Math.Cos(delta) * Math.Cos(alpha);
            LOS[1] = Math.Cos(delta) * Math.Sin(alpha);
            LOS[2] = Math.Sin(delta);
            
            return LOS;
        }

        public static double[] RaDec2AzEl(double Ra, double Dec, double lat, double lon, double jTime)
        {
            Ra = Astrolib.Deg2Rad(Ra);
            Dec = Astrolib.Deg2Rad(Dec);
            lat = Astrolib.Deg2Rad(lat);
            lat = Astrolib.Deg2Rad(lat);

            double LSTrad = LocalSiderialTime(jTime, lon);

            double LHA = LSTrad - Ra;

            double El = Math.Asin(Math.Sin(lat) * Math.Sin(Dec) + Math.Cos(lat) * Math.Cos(Dec) * Math.Cos(LHA));
            double Az = Math.Acos((Math.Sin(Dec) - Math.Sin(El) * Math.Sin(lat)) / (Math.Cos(El) * Math.Cos(lat)));

            double[] res = new double[2];
            res[0] = Astrolib.Rad2Deg(Az);
            res[1] = Astrolib.Rad2Deg(El);

            return res;
        }

        public static double GregorianToJD(int year, int month, int day, int hour, int minute, double second)       // gregorian date to julian
        {
            double JD = 367 * year + (int)((7 * (year + (int)((month + 9) / 12))) / 4) + (int)((275 * month) / 9) + day + 1721013.5 + (((second / 60) + minute) / 60 + hour) / 24;

            return JD;
        }

        public static double LocalSiderialTime(double jd, double lon) // local siderial time (angle)
        {
            // compute local siderial time
            double gmst = Astrolib.GetGMST(jd);
            
            double lonHrs = Astrolib.Rad2Deg(lon) / 15;
            double LSTHrs = gmst + lonHrs;
            
            double LSTdeg = LSTHrs * 15;							// [deg] local siderial time
            
            double LSTrad = Astrolib.Deg2Rad(LSTdeg);

            return LSTrad;
            
        }

        public static double GetGMST(double jd) 
        {
            // compute greenwhich mean siderial time 
            double gmst = 18.697374558 + 24.06570982441908 * (jd - 2451545);			
            gmst = gmst / 24;																// [hrs] Greenwhich mean siderial time
            
            return gmst;
        }

        public static double GetJulianFromUnix(int unixSecs)        // julian time from unix timestamp
        {
            return ( unixSecs / 86400.0 ) + 2440587.5;
        }

        public static Vector<double> SphereToEquatorial(Vector<double> meas, Vector<double> llaMat, double jd)
        {
            // compute right ascension and declination from azimuth and elevation measurements
            
            var V = Vector<double>.Build; 
            var RA_Dec = V.Dense(3);

            double az = meas[0];
            double el = meas[1];
            double lat = llaMat[0];
        
            double delta = Math.Asin(Math.Sin(el)*Math.Sin(lat) + Math.Cos(el)*Math.Cos(lat)*Math.Cos(az));         // declination
            double sinLHA = -((Math.Sin(az)*Math.Cos(el)) / Math.Cos(delta));
            double cosLHA = ((Math.Cos(lat)*Math.Sin(el) - Math.Sin(lat)*Math.Cos(az)*Math.Cos(el)) / Math.Cos(delta));
        
            double LHA = Math.Atan2(sinLHA, cosLHA);            // local hour angle
            double thetaLST = LocalSiderialTime(jd, llaMat[1]);         // local siderial time
            double alpha = thetaLST - LHA;                      // right ascension
        
            RA_Dec[0] = alpha;
            RA_Dec[1] = delta;
            
            return RA_Dec;
        }

        public static Vector<double> Cross(Vector<double> v1, Vector<double> v2)        // cross product
        {
            if ((v1.Count != 3 || v2.Count != 3))
            {
                string message = "Vectors must have a length of 3.";
                throw new Exception(message);
            }
            var Vx = Vector<double>.Build;
            var res = Vx.Dense(3);
            res[0] = v1[1] * v2[2] - v1[2] * v2[1];
            res[1] = -v1[0] * v2[2] + v1[2] * v2[0];
            res[2] = v1[0] * v2[1] - v1[1] * v2[0];

            return res;
        }

        public static Vector<double> Cartesian2Keplerian(Vector<double> rVec, Vector<double> vVec)  // cartesian coordinates to keplerian orbital elements
        {
            var V = Vector<double>.Build; 

            Vector<double> hVec = Astrolib.Cross(rVec, vVec);       // angular momentum vector

            double h = hVec.L2Norm();       // h magnitude
            double r = rVec.L2Norm();       // r magnitude
            double v = vVec.L2Norm();       // v magnitude

            Vector<double> kHat = V.Dense(3);       // kHat vector
            kHat[0] = 0;
            kHat[1] = 0;
            kHat[2] = 1;

            Vector<double> nVec = Astrolib.Cross(kHat, hVec);       // orbit normal vector
            // eccentricity vector (points to perigee)
            Vector<double> eVec = ((Math.Pow(v, 2) - mu / r) * rVec - (rVec.DotProduct(vVec) * vVec)) / mu;
            double e = eVec.L2Norm();       // eccentricity magnitude

            double eps = (Math.Pow(v, 2) - mu / r);     // orbital energy

            double sma;     // semi-major axis
            double p;       // semi-latus rectum
            if (e != 1.0) 
            {
                sma = -mu / (2 * eps);
                p = sma * (1 - Math.Pow(e, 2));
            } else {
                sma = double.PositiveInfinity;
                p = Math.Pow(h, 2) / mu;
            }

            double inc = Math.Acos(hVec[2] / h);        // orbit inclination (angle from equatorial plane)

            double node = Math.Acos(nVec[0] / nVec.L2Norm());       // right ascension of the ascending node (orientation of ascending node from first point of aries)
            if (nVec[1] < 0) {
                node = 2 * Math.PI - node;
            } 

            double argP = Math.Acos(nVec.DotProduct(eVec) / (nVec.L2Norm() * e));       // argument of perigee (orientation of perigee from line of nodes)
            if (eVec[2] < 0) {
                argP = 2 * Math.PI - argP;
            } 

            double truA = Math.Acos(eVec.DotProduct(rVec) / (e * r));       // true anomaly (angle of object with respect to perigee)
            if (rVec.DotProduct(vVec) < 0) {
                truA = 2 * Math.PI - truA;
            } 

            Vector<double> res = V.Dense(7);
            res[0] = sma;
            res[1] = e;
            res[2] = inc;
            res[3] = node;
            res[4] = argP;
            res[5] = truA;
            res[6] = p;

            return res;
        }

        public static Vector<double> GaussIOD(Matrix<double> Lmat, double[] jdVec, Matrix<double> rSiteMat) 
        {
            // Gauss method of angles only IOD
            // Uses three optical measurements to compute orbit
            
            double jd1 = jdVec[0];
            double jd2 = jdVec[1];
            double jd3 = jdVec[2];
            
            double tau1 = jd1 - jd2;        // time diff of 1 to 2
            double tau3 = jd3 - jd2;        // time diff of 2 to 3
            
            double a1 = tau3 / (tau3 - tau1);
            double a3 = -tau1 / (tau3 - tau1);
            double a1u = (tau3 * (Math.Pow((tau3 - tau1), 2) - Math.Pow(tau3, 2))) / (6 * (tau3 - tau1));
            double a3u = -(tau1 * (Math.Pow((tau3 - tau1), 2) - Math.Pow(tau1, 2))) / (6 * (tau3 - tau1));
            
            Matrix<double> Lmat_inv = Lmat.Inverse();
            Matrix<double> Mmat = Lmat_inv * rSiteMat;
            
            double d1 = Mmat[1, 0] * a1 - Mmat[1, 1] + Mmat[1, 2] * a3;
            double d2 = Mmat[1, 0] * a1u + Mmat[1, 2] * a3u;

            var V1 = Vector<double>.Build;
            var LmatSlice1 = V1.Dense(3);
            var LmatSlice2 = V1.Dense(3);
            var LmatSlice3 = V1.Dense(3);
            var rSiteMatSlice1 = V1.Dense(3);
            var rSiteMatSlice2 = V1.Dense(3);
            var rSiteMatSlice3 = V1.Dense(3);
            for (int i = 0; i < 3; i++)         // compute individual vectors (LOS, ground location) from matrices for calculations
            { 
                LmatSlice1[i] = Lmat[i, 0];
                LmatSlice2[i] = Lmat[i, 1];
                LmatSlice3[i] = Lmat[i, 2];
                rSiteMatSlice1[i] = rSiteMat[i, 0];
                rSiteMatSlice2[i] = rSiteMat[i, 1];
                rSiteMatSlice3[i] = rSiteMat[i, 2];
            }

            double C = LmatSlice1.DotProduct(rSiteMatSlice2);
            
            double rSite2_squared = Math.Pow(rSiteMatSlice2.L2Norm(), 2);
            
            double aPoly = -(Math.Pow(d1, 2) + 2 * C * d1 + rSite2_squared);        // first polynomial coefficient
            double bPoly = -(2 * mu * (C * d2 + d1 * d2));                          // second polynomial coefficient
            double cPoly = -Math.Pow(mu, 2) * Math.Pow(d2, 2);                      // third polynomial coefficient

            Polynomial p1 = new Polynomial(cPoly, 0, 0, bPoly, 0, 0, aPoly, 0, 1);
            Complex[] roots = p1.Roots();       // compute solutions to eighth order polynomial

            double bigr2 = -1e9;

            for (int i = 0; i < roots.Length; i++) {
                if (roots[i].IsReal() && roots[i].Re > bigr2) {
                    bigr2 = roots[i].Re;            // eliminate complex and nonpositive solutions (there should only be one positive real solution)
                }
            }
            
            double u = mu / Math.Pow(bigr2, 3);
            double c1 = a1 + a1u * u;
            double c2 = -1;
            double c3 = a3 + a3u * u;   


            var M1 = Matrix<double>.Build;
            var cVec = M1.Dense(3, 1);
            cVec[0, 0] = -c1;
            cVec[1, 0] = -c2;
            cVec[2, 0] = -c3;

            Matrix<double> rightSide = Mmat.Multiply(cVec);   

            var rhoVec = V1.Dense(3);               // slant ranges for each observation
            rhoVec[0] = rightSide[0, 0] / c1;
            rhoVec[1] = rightSide[1, 0] / c2;
            rhoVec[2] = rightSide[2, 0] / c3;

            var rhoVecLast = V1.Dense(3);
            for (int i = 0; i < 3; i++) {
                rhoVecLast[i] = rhoVec[i];
            }
            rhoVecLast[1] = -1e9;
            int iter = 0;

            Vector<double> r1;
            Vector<double> r2;
            Vector<double> r3;
            Vector<double> v2 = V1.Dense(3);
            while ((Math.Abs(rhoVecLast[1]-rhoVec[1]) > 1e-9))          // iterate over slant ranges until they stop changing
            {
                iter++;

                for (int i = 0; i < 3; i++) {
                    rhoVec[i] = rhoVecLast[i];
                }

                r1 = rhoVec[0] * LmatSlice1 + rSiteMatSlice1;
                r2 = rhoVec[1] * LmatSlice2 + rSiteMatSlice2;
                r3 = rhoVec[2] * LmatSlice3 + rSiteMatSlice3;

                v2 = Gibbs(r1, r2, r3, jd1, jd2, jd3);          // Gibbs method to find v2 using r1, r2, r3

                Vector<double> COE = Astrolib.Cartesian2Keplerian(r2, v2);      // get cartesian coordinates
                double p = COE[6];                              // semi-latus rectum

                double dAng1 = Astrolib.Angle1(r1, r2);         // angle between first and second positions
                double dAng3 = Astrolib.Angle1(r2, r3);         // angle betweeen second and third positions
                
                // f and g functions
                double f1 = 1 - (r1.L2Norm() / p) * (1 - Math.Cos(dAng1));  
                double f3 = 1 - (r3.L2Norm() / p) * (1 - Math.Cos(dAng3));
                double g1 = (r1.L2Norm() * r2.L2Norm() * Math.Asin(dAng1)) / Math.Sqrt(mu * p);
                double g3 = (r3.L2Norm() * r2.L2Norm() * Math.Asin(dAng3)) / Math.Sqrt(mu * p);

                c1 = g3 / (f1 * g3 - f3 * g1);
                c3 = -g1 / (f1 * g3 - f3 * g1);

                cVec[0, 0] = -c1;
                cVec[1, 0] = -c2;
                cVec[2, 0] = -c3;

                rightSide = Mmat.Multiply(cVec); 

                rhoVecLast[0] = rightSide[0, 0] / c1;           // recompute slant ranges
                rhoVecLast[1] = rightSide[1, 0] / c2;
                rhoVecLast[2] = rightSide[2, 0] / c3; 
            }
            
            r2 = rhoVec[1] * LmatSlice2 + rSiteMatSlice2;     // solution for r2

            Vector<double> stateVec = V1.Dense(6);
            for (int i = 0; i < 3; i++)
            {
                stateVec[i] = r2[i];
                stateVec[i + 3] = v2[i];
            }

            return stateVec;
            //Vector<double> cTimesRho = Mmat *      
            
        }

        public static Vector<double> Unit(Vector<double> vIn)       // compute unit vector
        {   
            var VOut = Vector<double>.Build;
            var res = VOut.Dense(3);

            double norm = vIn.L2Norm();
            for (int i = 0; i < 3; i++) {
                res[i] = vIn[i] / norm;
            }

            return res;
        }

        public static double Angle1(Vector<double> v1, Vector<double> v2)       // angle between two vectors
        {

            return Math.Acos(v1.DotProduct(v2) / (v1.L2Norm() * v2.L2Norm()));
        }

        public static Vector<double> Gibbs(Vector<double> r1, Vector<double> r2, Vector<double> r3, double jd1, double jd2, double jd3)    // Gibbs method
        {

            Vector<double> Z12 = Astrolib.Cross(r1, r2);  
            Vector<double> Z23 = Astrolib.Cross(r2, r3); 
            Vector<double> Z31 = Astrolib.Cross(r3, r1); 

            double alpha_cop = Math.Asin(Z23.DotProduct(r1) / (Z23.L2Norm() * r1.L2Norm()));

            if (alpha_cop > Astrolib.Deg2Rad(2))        // Check that observed positions are coplanar
            {
                // not coplanar, solution will likely be inacurrate
                string message = "Vectors are not coplanar.";
                throw new Exception(message);
            }

            Vector<double> N = r1.L2Norm() * Z23 + r2.L2Norm() * Z31 + r3.L2Norm() * Z12;
            Vector<double> D = Z12 + Z23 + Z31;

            if ((D.L2Norm() < 1e-8) || (N.L2Norm() < 1e-8) || (Astrolib.Unit(N).DotProduct(Astrolib.Unit(D)) < 1e-8))
            {
                return Astrolib.HerrickGibbs(r1, r2, r3, jd1, jd2, jd3);        // use Herrick Gibbs method if angles are small
            }

            Vector<double> S = (r2.L2Norm() - r3.L2Norm()) * r1 + (r3.L2Norm() - r1.L2Norm()) * r2 + (r1.L2Norm() - r2.L2Norm()) * r3;
            Vector<double> B = Astrolib.Cross(D, r2);

            double Lg = Math.Sqrt(mu / (N.L2Norm() * D.L2Norm()));
            
            Vector<double> v2 = (Lg / r2.L2Norm()) * B + Lg * S;        // v2 solution

            return v2;
        }
        public static Vector<double> HerrickGibbs(Vector<double> r1, Vector<double> r2, Vector<double> r3, double jd1, double jd2, double jd3)
        {
            // Herrick Gibbs method for small angles between observations

            double deltat31 = jd3 - jd1;
            double deltat32 = jd3 - jd2;
            double deltat21 = jd2 - jd1;

            Vector<double> v2 = -deltat32 * (1 / (deltat21 * deltat31) + mu / (12 * Math.Pow(r1.L2Norm(), 3))) * r1 + 
                                (deltat32 - deltat21) * (1 / (deltat21 * deltat32) + mu / (12 * Math.Pow(r2.L2Norm(), 3))) * r2 + 
                                deltat21 * (1 / (deltat32 * deltat31) + mu / (12 * Math.Pow(r3.L2Norm(), 3))) * r3;

            return v2;
        }   


    }

}