using MathNet.Numerics.LinearAlgebra;
using System;
using System.IO;
using System.Linq;

namespace AstrolibNS.IODNS
{
    public class IOD
    {

        static int[] unixTime = new int[3];         // unix timestamp
        static double[] latVec = new double[3];     // lat of observer
        static double[] lonVec = new double[3];     // lon of observer
        static double[] altVec = new double[3];     // alt of observer
        static double[] azVec = new double[3];      // azimuth measurement
        static double[] elVec = new double[3];      // elevation measurement
        public double[] stateVec = new double[6];   // output state vector

        public IOD(string dataFile)
        {
            IODdata(dataFile);          // Read observation data
            
            double[] jdVec = new double[3];
            var V = Vector<double>.Build;
            var M = Matrix<double>.Build;
            var llaMat = M.Dense(3, 3);             // latitude, longitude, and altitude of observations
            var measMat = M.Dense(2, 3);            // matrix of azimuth and elevation measurements
            for (int i = 0; i<3; i++) 
            {
                jdVec[i] = Astrolib.GetJulianFromUnix(unixTime[i]);     // vector of j2k seconds
                
                llaMat[0, i] = latVec[i];
                llaMat[1, i] = lonVec[i];
                llaMat[2, i] = altVec[i];
                
                measMat[0, i] = azVec[i];
                measMat[1, i] = elVec[i];
            }

            jdVec[0] = 2456159.986453;
            jdVec[1] = 2456159.991991;
            jdVec[2] = 2456159.994769;
            
            // loop three times
            var rSiteMat = M.Dense(3, 3);
            var LOS = M.Dense(3, 3);

            var rSite = V.Dense(3);
            var LOS_vec = V.Dense(3);

            for (int i = 0; i < 3; i++)
            {
                var V1 = Vector<double>.Build;
                var lla = V1.Dense(3);
                for (int j = 0; j < 3; j++)
                {
                    lla[j] = llaMat[j, i];
                }
                var meas = V1.Dense(2);
                for (int j = 0; j < 2; j++)
                {
                    meas[j] = measMat[j, i];
                }
                var jd = jdVec[i];

                rSite = Astrolib.GetGroundSiteMat(lla, jd);         // ECI coordinates of observation site
                Vector<double> RA_Dec = Astrolib.SphereToEquatorial(meas, lla, jd);     // right ascension and declination of observartions
                LOS_vec = Astrolib.GetLOSMat(RA_Dec);               // line of sight vector

                for (int j = 0; j < 3; j++)
                {
                    rSiteMat[j, i] = rSite[j];
                    LOS[j, i] = LOS_vec[j];
                }
            }

            Vector<double> res = Astrolib.GaussIOD(LOS, jdVec, rSiteMat);   // Gauss angles-only IOD method

            stateVec = res.ToArray();
        }

        private static void IODdata(string dataFile)
        {

            using (TextReader reader = File.OpenText("TestDataISS.txt"))
            {
                for (int i = 0; i < 3; i++) {
                    string text = reader.ReadLine();
                    string[] words = text.Split(' ');

                    unixTime[i] = int.Parse(words[1]);

                    latVec[i] = Astrolib.Deg2Rad(double.Parse(words[3]));
                    lonVec[i] = Astrolib.Deg2Rad(double.Parse(words[4]));
                    altVec[i] = Astrolib.Deg2Rad(double.Parse(words[5]));

                    azVec[i] = Astrolib.Deg2Rad(double.Parse(words[7]));
                    elVec[i] = Astrolib.Deg2Rad(double.Parse(words[9]));
                }
            }
        }
    }
}