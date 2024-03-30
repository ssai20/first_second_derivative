import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
public class CubicSplineError {
    public static void main(String[] args) {
        // задаём функцию и сетку
        for (int N = 16; N <= 256; N *= 2) {
            double[] h = new double[N]; // шаги
            double[] x = new double[N]; // точки сетки
            double[] y = new double[N]; // значения функции на сетке
            x[0] = 0.;
            for (int i = 0; i < N; i++) {
                h[i] = 1. / N;
                x[i+1] = x[i] + h[i]; // равномерная сетка
                y[i] = Math.sin(x[i]); // функция
            }

            // создаём объект кубического сплайна
            SplineInterpolator interpolator = new SplineInterpolator();
            PolynomialSplineFunction spline = interpolator.interpolate(x, y);

            // вычисляем погрешность в L точках
            int L = 5 * N; // количество точек, в которых будем вычислять погрешность
            double maxError = 0.0; // максимальная норма погрешности
            for (int i = 0; i < L; i++) {
                double xi = 0. + i * 1. / (L - 1); // точки, в которых вычисляем погрешность
                double error = Math.abs(spline.value(xi) - Math.sin(xi)); // погрешность
                if (error > maxError) {
                    maxError = error; // обновляем максимальную погрешность
                }
            }

            System.out.println("N = " + N + ", Максимальная норма погрешности: " + maxError);
        }
    }
}
