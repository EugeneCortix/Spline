﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using Microsoft.Win32;
using System.IO;

namespace Spline
{
    /// <summary>
    /// Логика взаимодействия для MainWindow.xaml
    /// </summary>

    

    public partial class MainWindow : Window
    {
        static int n;
        static List<double> x;
        static List<double> f;
        static List<double> P;
        double[] b;
        double[] q;
        double[,] A;
        List<double[,]> Aloc;
        List<double[]> bloc;
        List<double> xk;
        double xmin;
        double xmax;

        int al = 0;
        int li = 0;
        int ve = 0;

        public MainWindow()
        {
            InitializeComponent();
            
        }

        private void Button_Click(object sender, RoutedEventArgs e)
        {
            x = new List<double>();
            f = new List<double>();

            OpenFileDialog openFileDialog = new OpenFileDialog();
            openFileDialog.Title = "Open File";
            if (openFileDialog.ShowDialog() == false)
                return;
            string[] fileText = System.IO.File.ReadAllLines(openFileDialog.FileName);
            n = Convert.ToInt32(fileText[0]);
            for (int i = 1; i < n + 1; i++)
            {
                string[] xy;
                xy = fileText[i].Split('\t');
                f.Add(double.Parse(xy[1]));
                x.Add(double.Parse(xy[0]));
            }
            xmaxmins();
            assemblyA();
            assemblyb();
            q = solveMatrix();
            printvect(q, "q-vector");
            assemblyP();
            paintPoints();
            PrintGlob();
        }

        private void xmaxmins()
        {
            xmin = x[0];
            xmax = x[0];
            for (int i = 1; i < x.Count; i++)
            {
                if (xmin > x[i]) xmin = x[i];
                if (xmax < x[i]) xmax = x[i];
            }
        }

        private void paintPoints()
        {
            
            double ymin = f[0];
            double ymax = f[0];
            double ySubZero = 0;
            double xSubZero = 0;

            for (int i = 1; i < f.Count; i++)
            {
                if (ymin > f[i]) ymin = f[i];
                if (ymax < f[i]) ymax = f[i];
                if (ymin > P[i]) ymin = P[i];
                if (ymax < P[i]) ymax = P[i];
            }



            if (xmin < 0) xSubZero = xmin;
            if (ymin < 0) ySubZero = ymin;

            double kx = graphCanvas.ActualWidth /( xmax - xmin);
            double ky = graphCanvas.ActualHeight / (ymax - ymin);

            // Draw points

            for (int i = 0; i < x.Count; i++)
            {
                double px = (x[i] - xSubZero) * kx;
                double pfy = (f[i] - ySubZero -ymin) * ky;
                double ppy = (P[i] - ySubZero -ymin) * ky;

                
                Point point = new Point(px, 300 - pfy);
                Ellipse elipse = new Ellipse();

                elipse.Width = 4;
                elipse.Height = 4;

                elipse.StrokeThickness = 2;
                elipse.Stroke = Brushes.Red;
                elipse.Margin = new Thickness(point.X - 2, point.Y - 2, 0, 0);

                graphCanvas.Children.Add(elipse);
            }

            for (int i = 1; i < x.Count; i++)
            {
                double px = (x[i] - xSubZero) * kx;
                double ppy = (P[i] - ySubZero - ymin) * ky;
                double px2 = (x[i-1] - xSubZero) * kx;
                double ppy2 = (P[i-1] - ySubZero - ymin) * ky;



                Line line = new Line()
                {
                    X1 = px,
                    Y1 = 300 - ppy,
                    X2 = px2,
                    Y2 = 300 - ppy2,
                    Stroke = Brushes.Black,
                    StrokeThickness = 3
                }; 

                graphCanvas.Children.Add(line);
            }
        }
        private void drawSpline()
        {

        }
        private void assemblyA()
        {
            Aloc = new List<double[,]>();
            xk = new List<double>();
            // Create ranges
            xk.Add(xmin - Math.Abs(xmin * 0.2));
            // Crash into a parts
            int a = 2;
            for (int i = 1; i < a; i++)
            {
                xk.Add(xmin + (xmax - xmin)/a*i);
            }
            xk.Add(xmax + Math.Abs(xmax*0.2));
            printlist(xk, "xk");
            // Making local matrixes
            /*for(int i = 0; i < x.Count; i++)
            {
                int ind = getxk(x[i]);
                double eps = (x[i] - xk[ind]) / (xk[ind + 1] - xk[ind]);
                double[,] locmat = new double[4, 4];
                for (int k = 0; k < 4; k++)
                {
                    for (int l = 0; l < 4; l++)
                        locmat[k, l] = phi(eps, k + 1) * phi(eps, l + 1);
                }
                Aloc.Add(locmat);
            }*/

            for (int i = 1; i < xk.Count; i++)
            {
                // Create matrix
                double[,] locmat = new double[4, 4];
                // Find all x in this range
                List<double> xInArea = new List<double>();
                for(int s = 0; s< x.Count; s++)
                {
                    if (x[s] <= xk[i] && x[s] > xk[i - 1])
                        xInArea.Add(x[s]);
                }

                // Add elements
                for(int t = 0; t < xInArea.Count; t++)
                {
                   double eps = (xInArea[t] - xk[i - 1]) / (xk[i] - xk[i-1]);

                    for(int m = 0; m < 4; m++)
                        for(int l = 0; l < 4; l++)
                        {
                            locmat[m, l] += phi(eps, m + 1) * phi(eps, l + 1);
                        }
                }
                printmat(locmat, "Alocal");
                Aloc.Add(locmat);
            }

                // Global matrix
                int size = 2 + Aloc.Count*2;
            A = new double[size, size];
            for (int u = 0; u < Aloc.Count; u++)
            {
                for (int i = 0; i < 4; i++)
                {
                    for (int j = 0; j < 4; j++)
                    {
                        A[i + u*2, j + u*2] += Aloc[u][i,j];
                    }
                }
            }
            // PrintGlob();
            printmat(A, "Aglobal");
        }
        //Find xk
        private int getxk(double x)
        {
            int i = 0;
            while (xk[i] < x) i++;
            return i - 1 ;
        }
        // Printing the matrix and vectors to the text file
        private void PrintGlob()
        {
            int w = Convert.ToInt32(Math.Pow(A.Length, 0.5));
            string s = "Matrix A\n";
            for (int i = 0; i < w; i++)
            {
                for (int j = 0; j < 6; j++)
                {
                    s+=A[i, j].ToString()+'\t';
                }
                s += '\n';
            }
            s += "\n Vectors \n b\tq\n";
            for (int i = 0; i < w; i++)
            {
                s += b[i].ToString()+ '\t' + q[i].ToString() + '\n';
            }
            s += "P-vector\n";
            for (int i = 0; i < P.Count; i++)
            {
                s += P[i].ToString() + '\n';
            }
            System.IO.File.WriteAllText("..\\...\\...\\globals.txt", s);
        }
        // All phi-functions
        private double phi(double eps, int a)
        {
            switch (a)
            {
                case 1:
                    return 1 - 3 * Math.Pow(eps, 2) + 2 * Math.Pow(eps, 3);
                case 2:
                    return eps - 2 * Math.Pow(eps, 2) + Math.Pow(eps, 3);
                case 3:
                    return 3 * Math.Pow(eps, 2) - 2 * Math.Pow(eps, 3);
                case 4:
                    return -Math.Pow(eps, 2) +  Math.Pow(eps, 3);
            }
                
            return -1;
        }
        private void assemblyb()
        {
            // Making local vectors
            bloc = new List<double[]>();
            /*  for (int n = 0; n < x.Count; n++)
              {
                  int ind = getxk(x[n]);
                  double eps = (x[n] - xk[ind]) / (xk[ind + 1] - xk[ind]);
                  double[] locvect = new double[4];
                  for (int i = 0; i < 4; i++)
                  {
                      locvect[i] = f[n] * phi(eps, i + 1); 
                  }
                  bloc.Add(locvect);
              }*/
            for (int i = 1; i < xk.Count; i++)
            {
                List<double> xInArea = new List<double>();
                List<double> fInArea = new List<double>();
                double[] locvect = new double[4];
                for (int s = 0; s < x.Count; s++)
                {
                    if (x[s] <= xk[i] && x[s] > xk[i - 1])
                        xInArea.Add(x[s]);
                        fInArea.Add(f[s]);
                }

                // Add elements
                for (int t = 0; t < xInArea.Count; t++)
                {
                    double eps = (xInArea[t] - xk[i - 1]) / (xk[i] - xk[i - 1]);

                    for (int m = 0; m < 4; m++)
                        {
                        locvect[m] += phi(eps, m + 1) * fInArea[t];// fInArea
                        }
                }
                printvect(locvect,"blocal");
                bloc.Add(locvect);
            }

                // Creating of the global vector b
                b = new double[(xk.Count - 1)*2 +2]; 
            for (int n = 0; n < bloc.Count; n++)
            {
                for (int i = 0; i < 4; i++)
                    b[i + n * 2] = bloc[n][i];
            }
            printvect(b, "bglobal");
        }
        private void assemblyP()
        {
            P = new List<double>();
            for(int i = 0; i < x.Count; i++)
            {
                int ind = getxk(x[i]);
                double p = 0;
                double eps = (x[i] - xk[ind]) / (xk[ind + 1] - xk[ind]);
                for(int j = 0; j < 4; j++)
                {
                    p += q[ind * 2 + j] * phi(eps, j + 1);
                }
                P.Add(p);
            }
            printlist(P, "P-vector");
        }
        // Solving Aq=b equation
        private double[] solveMatrix()
        {
            double[] ans = new double[b.Length];
            solveSLAE(ans);
            return ans;
        }
        // Решает СЛАУ
        public void solveSLAE(double[] ans)
        {
            int nSLAE = b.Length;
            if (ans.Length != nSLAE)
                throw new Exception("Size of the input array is not compatable with size of SLAE");




            for (int i = 0; i < nSLAE; i++)
            {
                double del = A[i, i];
                double absDel = Math.Abs(del);
                int iSwap = i;


                for (int j = i + 1; j < nSLAE; j++) // ищем максимальный элемент по столбцу
                {
                    if (absDel < Math.Abs(A[j, i]))
                    {
                        del = A[j, i];
                        absDel = Math.Abs(del);
                        iSwap = j;
                    }
                }

                if (iSwap != i)
                {
                    double buf;
                    for (int j = i; j < nSLAE; j++)
                    {
                        buf = A[i, j];
                        A[i, j] = A[iSwap, j];
                        A[iSwap, j] = buf;
                    }
                    buf = b[i];
                    b[i] = b[iSwap];
                    b[iSwap] = buf;
                }

                for (int j = i; j < nSLAE; j++)
                    A[i, j] /= del;

                b[i] /= del;

                for (int j = i + 1; j < nSLAE; j++)
                {
                    if (A[j, i] == 0) continue;

                    double el = A[j, i];
                    for (int k = i; k < nSLAE; k++)
                    {
                        A[j, k] -= A[i, k] * el;
                    }

                    b[j] -= b[i] * el;
                }
            }

            for (int i = nSLAE - 1; i > -1; i--)
            {
                for (int j = i + 1; j < nSLAE; j++)
                    b[i] -= ans[j] * A[i, j];
                ans[i] = b[i];
            }
        }


        //Danila print 
        private void printvect(double[] vect, string vectname)
        {
            string s = "";
            for(int i = 0; i < vect.Length; i++)
            {
                s += vect[i].ToString() + '\n';
            }

            System.IO.File.WriteAllText("..\\...\\...\\" + vectname +ve.ToString()+ ".txt", s);
            ve++;
        }
        private void printmat(double[,] mat, string matname)
        {
            int w = Convert.ToInt32(Math.Pow(mat.Length, 0.5));
            string s = "";
            for (int i = 0; i < w; i++)
            {
                for (int j = 0; j < w; j++)
                {
                    s += mat[i, j].ToString() + '\t';
                }
                s += '\n';
            }
            System.IO.File.WriteAllText("..\\...\\...\\"+matname +al.ToString()+".txt", s);
            al++;
        }

        private void printlist(List<double> list, string lname)
        {
            string s = "";
            for (int i = 0; i < list.Count; i++)
            {
                s += list[i].ToString() + '\n';
            }

            System.IO.File.WriteAllText("..\\...\\...\\" + lname + li.ToString()+ ".txt", s);
            li++;
        }
    }
}
