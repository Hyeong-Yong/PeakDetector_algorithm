using MathNet.Numerics;
using MathNet.Numerics.Distributions;

namespace PeakDetector {
    public partial class Form1 : Form {
        public Form1() {
            InitializeComponent();
        }

        private void formsPlot1_Load(object sender, EventArgs e) {

        }

        SavitzkyGolayFilter sgFilter = new SavitzkyGolayFilter(3, 5);
        private void Form1_Load(object sender, EventArgs e) {
            double[] dataX = { 2800,    2801,   2802,   2803,   2804,   2805,   2806,   2807,   2808,   2809,   2810,   2811,   2812,   2813,   2814,   2815,   2816,   2817,   2818,   2819,   2820,   2821,   2822,   2823,   2824,   2825,   2826,   2827,   2828,   2829,   2830,   2831,   2832,   2833,   2834,   2835,   2836,   2837,   2838,   2839,   2840,   2841,   2842,   2843,   2844,   2845,   2846,   2847,   2848,   2849,   2850,   2851,   2852,   2853,   2854,   2855,   2856,   2857,   2858,   2859,   2860,   2861,   2862,   2863,   2864,   2865,   2866,   2867,   2868,   2869,   2870,   2871,   2872,   2873,   2874,   2875,   2876,   2877,   2878,   2879,   2880,   2881,   2882,   2883,   2884,   2885,   2886,   2887,   2888,   2889,   2890,   2891,   2892,   2893,   2894,   2895,   2896,   2897,   2898,   2899,   2900,   2901,   2902,   2903,   2904,   2905,   2906,   2907,   2908,   2909,   2910,   2911,   2912,   2913,   2914,   2915,   2916,   2917,   2918,   2919,   2920, };
            double[] dataY = {5.43703,  5.434557,   5.433024,   5.431118,   5.428057,    5.429489,   5.424885,   5.426195,   5.422085, 5.423,   5.427169,   5.431213,   5.429752,   5.431179,   5.433391,   5.437109,   5.435175,   5.43678,    5.435597,   5.436623,   5.435419,   5.435282,   5.435494,   5.434388,   5.436521,   5.437857,   5.435455,   5.435553,   5.430678,   5.434578,   5.432785,   5.43037,    5.428858,   5.421721,   5.422452,   5.422357,   5.421137,   5.41531,    5.407195,   5.402058,   5.397178,   5.389384,   5.379775,   5.375841,   5.37286,    5.36471,    5.354129,   5.336562,   5.324929,   5.316673,   5.300639,   5.286653,   5.265749,   5.245448,   5.213737,   5.18776,    5.146694,   5.111421,   5.057123,   5.026265,   5.005996,   4.997886,   4.996568,   4.999345,   4.995416,   4.971701,   4.987401,   5.004643,   5.036184,   5.069792,   5.105592,   5.140934,   5.133518,   5.192511,   5.228797,   5.253797,   5.264524,   5.276324,   5.295041,   5.301136,   5.32151,    5.341326,   5.348882,   5.356319,   5.360539,   5.367904,   5.374137,   5.38663,    5.393793,   5.400567,   5.407996,   5.413598,   5.412528,   5.413782,   5.423967,   5.427044,   5.42588,    5.427079,   5.426817,   5.428431,   5.431079,   5.430982,   5.430384,   5.429197,   5.430747,   5.428452,   5.427179,   5.423908,   5.425827,   5.422919,   5.426487,   5.425875,   5.425067,   5.423061,   5.424894,   5.426446,   5.42834,    5.425823,   5.42737,    5.428583,   5.431241};
            var data1= formsPlot1.Plot.Add.Scatter(dataX, dataY);
            data1.LegendText = "Origin";

            double[] dataY2 = PeakDect.MovingAverage(dataY, 3);
            //var data2 = formsPlot1.Plot.Add.Scatter(dataX, dataY2);
            //data2.LegendText = "Moving Avg";
            
            double[] dataY3 = sgFilter.Process(dataY);
            //var data3 = formsPlot1.Plot.Add.Scatter(dataX, dataY3);
            //data3.LegendText = "SavitzkyGolay";
            formsPlot1.Plot.ShowLegend();
            formsPlot1.Refresh();
            List<int> valley_data1 = PeakDect.MultipleValleyFinding(dataY);
            List<int> valley_data2 = PeakDect.MultipleValleyFinding(dataY2);
            List<int> valley_data3 = PeakDect.MultipleValleyFinding(dataY3);

            Console.WriteLine("Valley Detector: Origin");
            foreach (var value in valley_data1) {
                Console.WriteLine(dataX[value]);
            }
            Console.WriteLine("------");

            Console.WriteLine("Valley Detector: Moving Avg Filter ");
            foreach (var value in valley_data2) {
                Console.WriteLine(dataX[value]);
            }
            Console.WriteLine("------");

            Console.WriteLine("Valley Detector: Savitzky-Golay Filter ");
            foreach (var value in valley_data3) {
                Console.WriteLine(dataX[value]);
            }
            Console.WriteLine("------");


            double real_¥ò = 2.0;
            double real_¥ì = 2870;

            //Define gaussian function
            var gaussian = new Func<double, double, double, double>((¥ò, ¥ì, x) =>
            {
                return -Normal.PDF(¥ì, ¥ò, x)+5.45;
            });

            //Generate sample gaussian data
            var data = Enumerable.Range(0, 300)
                .Select(x => 2800 + x * 0.5)
                .Select(x => new { x, y = gaussian(real_¥ò, real_¥ì, x) });

            var xs = data.Select(d => d.x).ToArray();
            var ys = data.Select(d => d.y).ToArray();
            var initialGuess_¥ò = 2.0;
            var initialGuess_¥ì = 2870;

            var fit = Fit.Curve(dataX, dataY, gaussian, initialGuess_¥ò, initialGuess_¥ì);
            //fit.Item1 = ¥ò, fit.Item2 = ¥ì
            Console.WriteLine(fit.P0 +"//"+ fit.P1);
            var GaussianOrigin = formsPlot1.Plot.Add.Scatter(xs, ys);
            GaussianOrigin.LegendText = "Gauss Origin";
            //var Gaussianfit = formsPlot1.Plot.Add.Scatter(xs, fit);
            //Gaussianfit.LegendText = "Gauss Origin";
        }
    }
}
