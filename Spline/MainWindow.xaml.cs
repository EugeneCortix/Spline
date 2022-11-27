using System;
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
        static List<float> x;
        static List<float> y;
        double[] b = { 5, 6, 8, 9, 6 };
        double[,] A;
        List<double[,]> Aloc;
        float xmin;
        float xmax;

        public MainWindow()
        {
            InitializeComponent();
            solveMatrix();
        }

        private void Button_Click(object sender, RoutedEventArgs e)
        {
            x = new List<float>();
            y = new List<float>();

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
                y.Add(float.Parse(xy[1]));
                x.Add(float.Parse(xy[0]));
            }
            paintPoints();
            assemblyA();
        }

        private void paintPoints()
        {
            xmin = x[0];
            xmax = x[0];
            float xSubZero = 0;
            float ymin = y[0];
            float ymax = y[0];
            float ySubZero = 0;

            for (int i = 1; i < x.Count; i++)
            {
                if (xmin > x[i]) xmin = x[i];
                if (ymin > y[i]) ymin = y[i];
                if (xmax < x[i]) xmax = x[i];
                if (ymax < y[i]) ymax = y[i];
            }

            if (xmin < 0) xSubZero = xmin;
            if (ymin < 0) ySubZero = ymin;

            double kx = graphCanvas.ActualWidth /( xmax - xmin);
            double ky = graphCanvas.ActualHeight / (ymax - ymin);

            // Draw points

            for (int i = 0; i < x.Count; i++)
            {
                double px = (x[i] - xSubZero) * kx;
                double py = (y[i] - ySubZero) * ky;

                Point point = new Point(px, 300 - py);
                Ellipse elipse = new Ellipse();

                elipse.Width = 4;
                elipse.Height = 4;

                elipse.StrokeThickness = 2;
                elipse.Stroke = Brushes.Red;
                elipse.Margin = new Thickness(point.X - 2, point.Y - 2, 0, 0);

                graphCanvas.Children.Add(elipse);
            }
        }
        private void assemblyA()
        {
            Aloc = new List<double[,]>();
            List<double> xk = new List<double>();
            // Create ranges
            xk.Add(xmin - Math.Abs(xmin * 0.2));
            // Побить по-другому
            for (int i = 1; i < x.Count; i++)
            {
                xk.Add(x[i - 1] + (x[i] - x[i - 1]) / 2);
            }
            xk.Add(xmax + Math.Abs(xmax*0.2));

            // Making local matrixes
            for(int i = 0; i < x.Count; i++)
            {
                double eps = (x[i] - xk[i]) / (xk[i + 1] - xk[i]);
                double[,] locmat = new double[4, 4];
                for (int k = 0; k < 4; k++)
                {
                    for (int l = 0; l < 4; l++)
                        locmat[k, l] = phi(eps, k + 1) * phi(eps, l + 1);
                }
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
            PrintGlob();
        }
        private void PrintGlob()
        {
            string s = "";
            for (int i = 0; i < 6; i++)
            {
                for (int j = 0; j < 6; j++)
                {
                    s+=A[i, j].ToString()+'\t';
                }
                s += '\n';
            }
            System.IO.File.WriteAllText("..\\...\\...\\globalmatrix.txt", s);
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
                    return -Math.Pow(eps, 2) + 3 * Math.Pow(eps, 3);
            }
                
            return -1;
        }

        private void solveMatrix()
        {
            double[,] Matrix = { { 1, 2, 3, 4, 5 }, { 16, 27, 8,9,25 },{19,32,13,14,15},{16,17,19,89,20}, {21,22,23,33,25 } };
             
            
            int size = Convert.ToInt32(Math.Pow((Matrix.Length), 0.5));
            Matrix = triangle(Matrix);
            //Bm = triangle(Bm);
            double[] q = new double[size];
    
            int s = 1;
            //Solving SoLAE 
            for (int i = 0; i<size; i++)
            {
                double sum = 0;
                for (int j = 0; j < i; j++)
                {
                    sum += Matrix[size - 1 - i, j] * q[j];
                }
                q[i] = (b[size -1-i] - sum)/ Matrix[size - 1 - i, i];
            }

        }
        private double[,] triangle(double[,] Matrix)
        {
            int size = Convert.ToInt32(Math.Pow((Matrix.Length), 0.5));

            for (int l = 0; l < size-1; l++)
            {

                for (int i = 1; i +l < size; i++)
                {
                    double k = Matrix[i+l, size - 1 -l] / Matrix[l, size - 1 - l];
                    b[l + i] -= b[l] * k;
                    for (int j = 0; j < size; j++)
                    {
                        Matrix[i+l, j] -= Matrix[l, j] * k;
                    }
                }

            }

            // Computational error
            
            double eps = 0.0000000001;
            for (int i = 0; i < size; i++)
                for (int j = 0; j < size; j++)
                    if (Math.Abs(Matrix[i, j]) < eps) Matrix[i, j] = 0;
            return Matrix;
        }
    }
}
