using MathNet.Numerics;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;
using System.Numerics;

namespace PeakDetector {

    public sealed class SavitzkyGolayFilter {
        private readonly int sidePoints;

        private Matrix<double> coefficients;

        public SavitzkyGolayFilter(int sidePoints, int polynomialOrder) {
            this.sidePoints = sidePoints;
            Design(polynomialOrder);
        }

        /// <summary>
        /// Smoothes the input samples.
        /// </summary>
        /// <param name="samples"></param>
        /// <returns></returns>
        public double[] Process(double[] samples) {
            int length = samples.Length;
            double[] output = new double[length];
            int frameSize = (sidePoints << 1) + 1;
            double[] frame = new double[frameSize];

            Array.Copy(samples, frame, frameSize);

            for (int i = 0; i < sidePoints; ++i) {
                output[i] = coefficients.Column(i).DotProduct(MathNet.Numerics.LinearAlgebra.Vector<double>.Build.DenseOfArray(frame));
            }

            for (int n = sidePoints; n < length - sidePoints; ++n) {
                Array.ConstrainedCopy(samples, n - sidePoints, frame, 0, frameSize);
                output[n] = coefficients.Column(sidePoints).DotProduct(MathNet.Numerics.LinearAlgebra.Vector<double>.Build.DenseOfArray(frame));
            }

            Array.ConstrainedCopy(samples, length - frameSize, frame, 0, frameSize);

            for (int i = 0; i < sidePoints; ++i) {
                output[length - sidePoints + i] = coefficients.Column(sidePoints + 1 + i).DotProduct(MathNet.Numerics.LinearAlgebra.Vector<double>.Build.Dense(frame));
            }

            return output;
        }

        private void Design(int polynomialOrder) {
            double[,] a = new double[(sidePoints << 1) + 1, polynomialOrder + 1];

            for (int m = -sidePoints; m <= sidePoints; ++m) {
                for (int i = 0; i <= polynomialOrder; ++i) {
                    a[m + sidePoints, i] = Math.Pow(m, i);
                }
            }

            Matrix<double> s = Matrix<double>.Build.DenseOfArray(a);
            coefficients = s.Multiply(s.TransposeThisAndMultiply(s).Inverse()).Multiply(s.Transpose());
        }
    }


    public class PeakDect{









        public static double[] MovingAverage(double[] source, int period) {
            var ma = new double[source.Length];
            double sum = 0;
            for (int bar = 0; bar < period; bar++)
                sum += source[bar];

            ma[period - 1] = sum / period;

            for (int bar = period; bar < source.Length; bar++)
                ma[bar] = ma[bar - 1] + source[bar] / period
                                      - source[bar - period] / period;
            
            for (int i = 0; i < period; i++) {
                ma[i] = source[i];
            }
            return ma;
        }
        public static double[] MovAvg_Filter(double[] sample, int nFilterSize) {
            List<double> lst = new List<double>(); // 결과로 반환할 데이터 배열 
            double[] arrNumbers = new double[3]; // 필터 배열 
            int pos = 0; //필터시작위치 
            double newAvg = 0; // 평균 
            double sum = 0; // 합계 
            int len = arrNumbers.Length; // 필터 크기 
            int count = sample.Length; // 입력받은 데이터 크기 

            for (int i = 0; i < count; i++) // 입력받은 데이터 만큼 반복 
            {
                sum = sum - arrNumbers[pos] + sample[i]; // 필터 내의 데이터와 현재 데이터의 합계를 구한다.
                arrNumbers[pos] = sample[i]; // 필터에 추가 
                newAvg = sum / len; // 평균 구하기 pos++;

                // 필터 증가 
                if (pos >= len) {
                    pos = 0; // 필터 처음부터 
                }

                // 필터가 커지면서 데이터가 뒤로 밀리는 것을 방지하기 위해 
                // 필터링 이전 값은 최초 필터링 된 평균값으로 대체합니다.

                if (i == nFilterSize) {
                    for (int j = 0; j < lst.Count; j++) {
                        lst[j] = newAvg;
                    }
                }

                // 결과에 추가
                lst.Add(newAvg);
            }

            return lst.ToArray(); // 결과 반환 
        }

        /// <summary>
        /// 누적 평균값을 구함(?) 이상함
        /// </summary>
        /// <param name="signal"></param>
        /// <returns></returns>
        public static List<double> AvgFilter(double[] signal) {
            double prev_average = 0;
            int num_of_data = 1;
            List<double> ret = new List<double>();
            foreach (var value in signal) {
                double average, alpha;
                //
                //alpha of average filter
                alpha = (num_of_data - 1) / (num_of_data + 0.0);
                //recursive expression of average filter
                average = alpha * prev_average + (1 - alpha) * value;
                //update previous state value of average filter, and +1 the number of sample
                prev_average = average;
                num_of_data += 1;
                ret.Add(average);
            }

            return ret;

        }

        /// <summary>
        /// <para>Implements a Savitzky-Golay smoothing filter, as found in [1].</para>
        /// <para>[1] Sophocles J.Orfanidis. 1995. Introduction to Signal Processing. Prentice-Hall, Inc., Upper Saddle River, NJ, USA.</para>
        /// </summary>




     public static int NaivePeakFinding(double[] signal) {

            int? peakIndex = null;
            double? peakValue = null;
            int index = 0;
            foreach (var value in signal) {
                if (peakValue == null || value > peakValue) {
                    peakIndex = index;
                    peakValue = value;
                }
                index++;
            }

            if (peakIndex == null) {
                return 0;
            }

            return (int)peakIndex;
        }
        public static int NaiveValleyFinding(double[] signal) {

            int? peakIndex = null;
            double? peakValue = null;
            int index = 0;
            foreach (var value in signal) {
                if (peakValue == null || value < peakValue) {
                    peakIndex = index;
                    peakValue = value;
                }
                index++;
            }

            if (peakIndex == null) {
                return 0;
            }

            return (int)peakIndex;
        }

        public static List<int> MultiplePeakFinding(double[] signal, bool isAvgFilter = true) {
            List<int> peakIndices = new List<int> { };
            double baseline;

            if (isAvgFilter == true) {
                baseline = AvgFilter(signal).Average();
            }
            else {
                baseline = signal.Average();
            }

            int? peakIndex = null;
            double? peakValue = null;
            int index = 0;

            foreach (var value in signal) {
                if (value > baseline) {
                    if (peakValue == null || value > peakValue) {
                        peakIndex = index;
                        peakValue = value;
                    }
                }
                else if (value < baseline && peakIndex != null) {
                    peakIndices.Add((int)peakIndex);
                    peakIndex = null;
                    peakValue = null;
                }
                index++;
            }
            if (peakIndex != null) {
                peakIndices.Add((int)peakIndex);
            }
            return peakIndices;
        }
        public static List<int> MultipleValleyFinding(double[] signal, bool isAvgFilter = true) {
            List<int> peakIndices = new List<int> { };
            double baseline;

            if (isAvgFilter == true) {
                baseline = AvgFilter(signal).Average();
            }
            else {
                baseline = signal.Average();
            }

            int? peakIndex = null;
            double? peakValue = null;
            int index = 0;

            foreach (var value in signal) {
                if (value < baseline) {
                    if (peakValue == null || value < peakValue) {
                        peakIndex = index;
                        peakValue = value;
                    }
                }
                else if (value > baseline && peakIndex != null) {
                    peakIndices.Add((int)peakIndex);
                    peakIndex = null;
                    peakValue = null;
                }
                index++;
            }
            if (peakIndex != null) {
                peakIndices.Add((int)peakIndex);
            }
            return peakIndices;
        }

    }

}

//algorithm NaivePeakFinding(signal):
//// INPUT
////    signal = an array of signal values
//// OUTPUT
////    index = the index of the peak value in the signal
//
//peakIndex < -null
//    peakValue < -null
//    for index, value in signal: 
//        if peakValue = null or value > peakValue:
//            peakIndex < -index
//            peakValue < -value
//    return peakIndex

//algorithm MultiplePeakFinding(signal):
//    // INPUT
//    //    signal = an array of signal values
//    // OUTPUT
//    //    indices = an array of indices of the peak values in the signal

//    peakIndices < -an empty array
//    peakIndex <- null
//    peakValue <- null
//    baseline <- average(signal) 
//    for index, value in signal: 
//        if value > baseline:
//            if peakValue = null or value > peakValue:
//                peakIndex < -index
//                peakValue < -value
//        else if value < baseline and peakIndex != null:
//            peakIndices.append(peakIndex)
//            peakIndex < -null
//            peakValue < -null
//    if peakIndex != null:
//        peakIndices.append(peakIndex)
//    return peakIndices

//algorithm NoisyPeakFinding(signal):
//    // INPUT
//    //   signal = an array of signal values
//    // OUTPUT
//    //   indices = an array of indices of the peak values in the signal

//    smoothed < - []
//    for index in signal:
//        smoothed.append(average(signal[index - 2 : index + 2]))


//    peakIndices < -an empty array
//    peakIndex < -null
//    peakValue < -null
//    baseline < -average(smoothed)
//    for index, value in smoothed: 
//        if value > baseline:
//            if peakValue = null or value > peakValue:
//                peakIndex < -index
//                peakValue < -value
//        else if value < baseline and peakIndex != null:
//            peakIndices.append(peakIndex)
//            peakIndex < -null
//            peakValue < -null
//    if peakIndex != null:
//        peakIndices.append(peakIndex)
//    return peakIndices


