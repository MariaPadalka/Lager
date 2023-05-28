using System;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Text;

namespace Laguerre
{
    public static class MathUtils
    {
        public static double[] Linspace(double a, double b, int number_of_points)
        {
            if (number_of_points < 2)
            {
                throw new ArgumentException("The number of points must be at least 2.", nameof(number_of_points));
            }

            double[] points = Enumerable.Range(0, number_of_points)
                .Select(i => a + i * (b - a) / (number_of_points - 1))
                .ToArray();

            return points;
        }
    }

    public static class IntegrationUtils
    {
        public static double[] ArrIntegrals(Func<double, double> f, double[] points, double delta, LaguerreTransform laguerreObject)
        {
            double[] result = Enumerable.Range(0, laguerreObject.N + 1)
                                         .Select(k =>
                                             points.Select(p => f(p) * laguerreObject.LaguerreFunc(k, p) * Math.Exp(-laguerreObject.Alpha * p) * delta)
                                                   .Sum()
                                         )
                                         .ToArray();
            return result;
        }

        public static double[] ArrIntegrals(Func<double, double, double, double> f, double mu, double lmbd, double[] points, double delta, LaguerreTransform laguerreObject)
        {
            double[] result = Enumerable.Range(0, laguerreObject.N + 1)
                                         .Select(k =>
                                             points.Select(p => f(p, mu, lmbd) * laguerreObject.LaguerreFunc(k, p) * Math.Exp(-laguerreObject.Alpha * p) * delta)
                                                   .Sum()
                                         )
                                         .ToArray();
            return result;
        }
    }

    public static class Function
    {
        static public double[] Func(double[] t)
        {
            double[] result = new double[t.Length];
            for (int i = 0; i < t.Length; i++)
            {
                double val = t[i];
                if (val >= 0 && val <= 2 * Math.PI)
                {
                    result[i] = Math.Sin(val - Math.PI / 2) + 1;
                }
                else if (val >= 2 * Math.PI)
                {
                    result[i] = 0;
                }
            }
            return result;
        }

        static public double Func(double t)
        {
            if (t >= 0 && t <= 2 * Math.PI)
            {
                return Math.Sin(t - Math.PI / 2) + 1;
            }
            else if (t >= 2 * Math.PI)
            {
                return 0;
            }
            return 0;
        }
    }

    public static class GaussianFunction
    {
        static public double[] Func(double[] t, double mu, double lmbd)
        {
            double[] result = new double[t.Length];
            for (int i = 0; i < t.Length; i++)
            {
                result[i] = (1 / (lmbd * Math.Sqrt(2 * Math.PI))) * Math.Exp(-(Math.Pow((t[i] - mu), 2)) / (2 * Math.Pow(lmbd, 2)));
            }
            return result;

        }
        static public double Func(double t, double mu, double lmbd)
        {
            var result = (1 / (lmbd * Math.Sqrt(2 * Math.PI))) * Math.Exp(-(Math.Pow((t - mu), 2)) / (2 * Math.Pow(lmbd, 2)));
            return result;
        }
    }


    public class LaguerreTransform
    {
        public LaguerreTransform(double beta, double sigma, double T_, int number_of_points, double eps, int N_)
        {
            Beta = beta;
            Sigma = sigma;
            Alpha = sigma - beta;
            NumberOfPoints = number_of_points;
            T = T_;
            Eps = eps;
            N = N_;
            InputArray = input_array();
        }

        public double Beta
        { get; set; }

        public double Sigma
        { get; set; }

        public double Alpha
        { get; set; }

        public double T
        { get; set; }

        public int NumberOfPoints
        { get; set; }

        public double Eps
        { get; set; }

        public int N
        { get; set; }

        public double[] InputArray
        { get; set; }


        public double[] input_array()
        {
            double[] input_array = new double[NumberOfPoints];
            double delta = T / (NumberOfPoints - 1);
            for (int i = 0; i < NumberOfPoints; i++)
            {
                input_array[i] = i * delta;
            }
            return input_array;
        }

        public double LaguerreFunc(int n, double t)
        {
            if (t < 0 && n < 0)
            {
                Console.WriteLine("You entered incorrect data.");
                return 0;
            }

            double first = Math.Sqrt(Sigma) * Math.Exp(t * (-Beta / 2));
            double second = first * (1 - (t * Sigma));

            if (n == 0)
            {
                return first;
            }
            else if (n == 1)
            {
                return second;
            }
            else
            {
                double third = 0;
                for (int i = 2; i <= n; i++)
                {
                    third = ((2 * i - 1 - t * Sigma) / i) * second - (i - 1) * first / i;
                    first = second;
                    second = third;
                }
                return third;
            }
        }

        public Tuple<double[], double[]> TabulationOfLaguerreFunction(int n)
        {
            double[] result_t = new double[NumberOfPoints];
            double[] lag = new double[NumberOfPoints];
            for (int i = 0; i < NumberOfPoints; i++)
            {
                result_t[i] = InputArray[i];
                lag[i] = LaguerreFunc(n, result_t[i]);
            }
            return Tuple.Create(result_t, lag);
        }

        public double[] LaguerreTransformation(Func<double, double> func)
        {
            int numberOfPoints = NumberOfPoints;
            double delta = (T - 0) / (numberOfPoints - 1);
            double halfDelta = delta / 2;
            double[] points = new double[numberOfPoints - 1];

            for (int i = 0; i < numberOfPoints - 1; i++)
            {
                points[i] = 0 + halfDelta + i * delta;
            }

            double[] funcRes0 = new double[N + 1];
            double[] funcRes1 = IntegrationUtils.ArrIntegrals(func, points, delta, this);

            while (funcRes0.Select((x, i) => Math.Abs(x - funcRes1[i])).Max() > Eps)
            {
                funcRes0 = funcRes1;
                numberOfPoints *= 2;
                double delta_ = (T - 0) / (numberOfPoints - 1);
                double halfDelta_ = delta_ / 2;
                double[] points_ = new double[numberOfPoints - 1];

                for (int i = 0; i < numberOfPoints - 1; i++)
                {
                    points_[i] = 0 + halfDelta_ + i * delta_;
                }

                funcRes1 = IntegrationUtils.ArrIntegrals(func, points, delta, this);
            }
            return funcRes1;
        }

        // added
        public double[] LaguerreTransformation(Func<double, double, double, double> func, double mu, double lmbd)
        {
            int numberOfPoints = NumberOfPoints;
            double delta = (T - 0) / (numberOfPoints - 1);
            double halfDelta = delta / 2;
            double[] points = new double[numberOfPoints - 1];

            for (int i = 0; i < numberOfPoints - 1; i++)
            {
                points[i] = 0 + halfDelta + i * delta;
            }

            double[] funcRes0 = new double[N + 1];
            double[] funcRes1 = IntegrationUtils.ArrIntegrals(func, mu, lmbd, points, delta, this);

            while (funcRes0.Select((x, i) => Math.Abs(x - funcRes1[i])).Max() > Eps)
            {
                funcRes0 = funcRes1;
                numberOfPoints *= 2;
                double delta_ = (T - 0) / (numberOfPoints - 1);
                double halfDelta_ = delta_ / 2;
                double[] points_ = new double[numberOfPoints - 1];

                for (int i = 0; i < numberOfPoints - 1; i++)
                {
                    points_[i] = 0 + halfDelta_ + i * delta_;
                }

                funcRes1 = IntegrationUtils.ArrIntegrals(func, mu, lmbd, points, delta, this);
            }
            return funcRes1;
        }


        public double InverseLaguerreTransform(double t, Func<double, double> func)
        {
            double[] sequence = LaguerreTransformation(func);
            double h = 0;
            for (int i = 0; i < sequence.Length; i++)
            {
                double lag = LaguerreFunc(i, t);
                h += sequence[i] * lag;
            }
            return h;
        }

        public double InverseLaguerreTransform(double t, Func<double, double, double, double> func, double mu, double lmbd)
        {
            double[] sequence = LaguerreTransformation(func, mu, lmbd);
            double h = 0;
            for (int i = 0; i < sequence.Length; i++)
            {
                double lag = LaguerreFunc(i, t);
                h += sequence[i] * lag;
            }
            return h;
        }
    }


    class Data
    {
        public static List<string> readDataFromFile(string path)
        {
            List<string> data = new List<string>();
            using (var file = new StreamReader(path))
            {
                string line;
                while ((line = file.ReadLine()) != null)
                {
                    data.Add(line);
                }
            }
            return data;
        }

        public static void writeTabulationDataToFile(string path, double res, List<Tuple<double[], double[]>> tabulation_data, Tuple<double[], double[]> sin_tab, double[] inv_x, double[] inv_gaus_x, List<List<double>> inverse, Tuple<double[], double[]> gaus_tab, List<List<double>> gaus_inverse, List<int> N_arr, Stopwatch stopwatch)
        {


            if (File.Exists(path)) File.Delete(path);
            // create the file

            var fs = new FileStream(path, FileMode.CreateNew);

            using (var writer = new StreamWriter(fs))
            {
                writer.Write("Calculate:");
                writer.Write($"\n{res}");


                writer.Write("\nTabulation:");
                writer.Write("\nx values\n");
                for (int i = 0; i < tabulation_data[0].Item1.Length; i++)
                {
                    writer.Write($"{tabulation_data[0].Item1[i]}\n");
                }


                writer.Write("y values\n");
                var count = 0;
                foreach (var t in tabulation_data)
                {
                    writer.Write($"n = {count}\n");
                    for (int i = 0; i < t.Item1.Length; i++)
                    {
                        writer.Write($"{t.Item2[i]}\n");
                    }
                    count++;
                }


                writer.Write("Sin\n");
                writer.Write("x\n");
                for (int i = 0; i < sin_tab.Item1.Length; i++)
                {
                    writer.Write($"{sin_tab.Item1[i]}\n");
                }

                writer.Write("y\n");
                for (int i = 0; i < sin_tab.Item2.Length; i++)
                {
                    writer.Write($"{sin_tab.Item2[i]}\n");
                }



                writer.Write("Gaussian\n");
                writer.Write("x\n");
                for (int i = 0; i < gaus_tab.Item1.Length; i++)
                {
                    writer.Write($"{gaus_tab.Item1[i]}\n");
                }

                writer.Write("y\n");
                for (int i = 0; i < gaus_tab.Item2.Length; i++)
                {
                    writer.Write($"{gaus_tab.Item2[i]}\n");
                }




                writer.Write("Inverse transformation:\n");
                writer.Write("x sin inverse\n");
                for (int i = 0; i < inv_x.Length; i++)
                {
                    writer.Write($"{inv_x[i]}\n");
                }


                writer.Write("y sin inverse\n");
                var c = 0;
                foreach (var inner in inverse)
                {
                    writer.Write($"n = {N_arr[c]}\n");
                    foreach (var x in inner)
                    {
                        writer.Write($"{x}\n");
                    }
                    c++;
                }


                ////// Gaussian
                writer.Write("x gaussian inverse\n");
                for (int i = 0; i < inv_gaus_x.Length; i++)
                {
                    writer.Write($"{inv_gaus_x[i]}\n");
                }


                writer.Write("y gaussian inverse\n");
                var o = 0;
                foreach (var inner in gaus_inverse)
                {
                    writer.Write($"n = {N_arr[o]}\n");
                    foreach (var x in inner)
                    {
                        writer.Write($"{x}\n");
                    }
                    o++;
                }


                writer.Write("Total running time:\n");
                writer.Write(stopwatch.Elapsed);

            }
        }
    }

    class MainClass
    {
        static public void Main(String[] args)
        {
            Stopwatch stopwatch = new Stopwatch();
            stopwatch.Start();

            NumberFormatInfo provider = new NumberFormatInfo();
            provider.NumberDecimalSeparator = ".";

            var data = Data.readDataFromFile("widgets_data.txt");

            double t = Convert.ToDouble(data[0], provider);
            double beta = Convert.ToDouble(data[1], provider);
            double sigma = Convert.ToDouble(data[2], provider);
            int N = int.Parse(data[3]);
            int n = int.Parse(data[4]);
            double T = Convert.ToDouble(data[5], provider);
            int number_of_points = int.Parse(data[6]);
            double eps = Convert.ToDouble(data[7], provider);
            double a = Convert.ToDouble(data[8], provider);
            double b = Convert.ToDouble(data[9], provider);
            int number_of_p = int.Parse(data[10]);
            

            string[] N_arr_str = data[11].Split(", ");
            List<int> N_arr = new List<int>();

            foreach (var k in N_arr_str)
            {
                N_arr.Add(int.Parse(k));
            }

            double a_ = Convert.ToDouble(data[12], provider);
            double b_ = Convert.ToDouble(data[13], provider);
            double mu = Convert.ToDouble(data[15], provider);
            double lambda = Convert.ToDouble(data[14], provider);



            var lag = new LaguerreTransform(beta, sigma, T, number_of_points, eps, N);


            // calculation
            var res = lag.LaguerreFunc(n, t);


            // tabulation
            List<Tuple<double[], double[]>> tabulation_data = new List<Tuple<double[], double[]>>();

            var i = 0;
            while (i < n)
            {
                var tabulation = lag.TabulationOfLaguerreFunction(i);

                tabulation_data.Add(tabulation);

                i++;
            }


            // tabulation sin [0; 2pi]
            var sequence = MathUtils.Linspace(a, b, number_of_p);
            var sin_res = Function.Func(sequence);
            var sin_tab = Tuple.Create(sequence, sin_res);

            // tabulation Gaussian function
            var sequence_g = MathUtils.Linspace(a_, b_, number_of_p);
            var gaus_res = GaussianFunction.Func(sequence_g, mu, lambda);
            var gaus_tab = Tuple.Create(sequence_g, gaus_res);


            // inverse
            List<List<double>> inverse = new List<List<double>>();
            List<List<double>> gaus_inverse = new List<List<double>>();


            var half_delta = ((b - a) / (number_of_p - 1) / 2);
            var t_range = MathUtils.Linspace(a + half_delta, b - half_delta, number_of_p - 1);
            var t_range2 = MathUtils.Linspace(a_ + half_delta, b_ - half_delta, number_of_p - 1);


            foreach (var y in N_arr)
            {
                lag.N = y;
                List<double> temp = new List<double>();
                List<double> gaus_temp = new List<double>();

                foreach (var g in t_range)
                {
                    temp.Add(lag.InverseLaguerreTransform(g, Function.Func));
                }
                foreach(var g in t_range2)
                {
                    gaus_temp.Add(lag.InverseLaguerreTransform(g, GaussianFunction.Func, mu, lambda));
                }
                inverse.Add(temp);
                gaus_inverse.Add(gaus_temp);
            }





            stopwatch.Stop();

            Data.writeTabulationDataToFile("calculatedData.txt", res, tabulation_data, sin_tab, t_range, t_range2, inverse, gaus_tab, gaus_inverse, N_arr, stopwatch);

            // Read the file contents into a string
            string fileContents = System.IO.File.ReadAllText("calculatedData.txt");
            // Replace all commas with dots in the string
            string modifiedContents = fileContents.Replace(",", ".");
            // Write the modified string back to the file
            System.IO.File.WriteAllText("calculatedData.txt", modifiedContents);

            Console.WriteLine("OLL IS WELL");
        }
    }
}