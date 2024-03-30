import java.util.*;
import java.util.function.BiFunction;
import java.util.function.Function;
//import org.apache.commons.math3.analysis.interpolation.NaturalSplineInterpolator;
import org.apache.commons.math3.analysis.interpolation.AkimaSplineInterpolator;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

public class Zadorin_1 {

    public static void findU (int N, double epsilon) {
        int L = 5*N;
        double[] U = new double[N+1];
        double[] proizvU = new double[L+1];
        double[] proizvPogranSloiU = new double[L+1];
        double h = 1/N;
        double hh = 1/L;
        double[] uzelX = new double[N+1];
        double[] x = new double[L+1];

        uzelX[0] = 0;
        for (int i=0;i<N;i++){
            uzelX[i] = uzelX[i-1]+h;
        }
        x[0] = 0;
        for (int j=1;j<L+1;j++){
            x[j] = x[j-1]+hh;
        }

        for (int i=0;i<N;i++){
            U[i] = Math.cos(Math.PI*uzelX[i]) + Math.exp(-Math.exp(uzelX[i]/epsilon));
        }
        for (int i=0;i<N;i++){
            proizvU[i] = (U[i] - U[i-1])/h;
            proizvPogranSloiU[i] = (-1/epsilon)*Math.exp(-uzelX[i]/epsilon)*(U[i] - U[i-1])/(Math.exp(-uzelX[i]/epsilon) - Math.exp(-uzelX[i-1]/epsilon));
        }
    }
    public static void progonka(double[] A, double[] B, double[] C, double[] F, int N, double[] c, double[] f, double[] h, double epsilon, int iteration){
        double[]alpha = new double[N+1];
        double[]betta = new double[N+1];


        double proizv1 = -1.0/epsilon; //производная Math.cos(Math.PI*x/2.) + Math.exp(-x/epsilon) в 0
        double proizv2 = -Math.PI/2.- Math.exp(-1.0/epsilon)/epsilon; //производная Math.cos(Math.PI*x/2.) + Math.exp(-x/epsilon) в 1

//        alpha[1] = h[1]/6.0;
//        betta[1] = proizv1 - (f[1]-f[0])/h[1];
        alpha[1] = 0.0;//мой метод
        betta[1] = -Math.PI*Math.PI/4.+1./(epsilon*epsilon);//мой метод
        for(int i=1;i<N;i++)
        {
            alpha[i+1]=B[i]/(C[i]-alpha[i]*A[i]);

            betta[i+1]=(A[i]*betta[i]+F[i])/(C[i]-alpha[i]*A[i]);
        }

//        c[0]=-Math.PI*Math.PI/4.0+1.0/(epsilon*epsilon); //вторая производная cos(Pi*x/4)+epsilon в 0
//        c[N]=Math.exp(-1./epsilon)/(epsilon*epsilon);//вторая производная cos(Pi*x/4)+epsilon в 1
//        c[0] = 0.;//sin^3
//        c[N] = 3.*(2.*Math.sin(1.)*Math.cos(1.)*Math.cos(1.) - Math.sin(1.)*Math.sin(1.)*Math.sin(1.));//sin^3
//        System.out.println(3.*(2.*Math.sin(1.)*Math.cos(1.)*Math.cos(1.) - Math.sin(1.)*Math.sin(1.)*Math.sin(1.)));




        c[0] = 6.*alpha[1]/(h[1]*(2.*alpha[1]+1.)) * ((h[1]*betta[1])/(6.*alpha[1]) + (f[1]-f[0])/h[1] - proizv1);
        c[N] = 6.*((f[N-1]-f[N])/h[N] - h[N]*betta[N]/6. + proizv2)/(h[N]*(2.+alpha[N]));


        for (int i=N-1;i>0;i--)
        {
            c[i]=alpha[i+1]*c[i+1] + betta[i+1];

            iteration++;
        }
    }

    public static  double[] ravnomSetka(int N){
//        System.out.println("Uniform mesh");
        double[] h = new double[N+1];
        for (int i=1;i<N+1;i++){
            h[i] = 1./N;
        }
        return h;
    }

    public static double[] setkaShishkina(double epsilon, int N){
//      System.out.println("Shishkin mesh");
        double sigma1;
        double[] h = new double[N+1];
        sigma1 = Math.min(0.5, 4.*epsilon*Math.log(N*1.0)/1.);
//            sigma2 = Math.min(0.5, 2.*epsilon*Math.log(N*1.0)/2.);

        for(int i=1;i<=N;i++){
            if ((i<=N/2)&&(i>=1)) {
                h[i] = 2.*sigma1/N;
//                    hk[i] = 2.*sigma2/N;
            }
            if ((i<=N)&&(i>N/2)) {
                h[i] = 2.*(1.-sigma1)/N;
            }
        }

        double uzelx[] = new double[N+1];
        uzelx[0]=0.;
//            uzely[0]=0.;
        for (int i=1;i<N+1;i++)
        {
            uzelx[i]=uzelx[i-1]+h[i];

        }
        uzelx[N/2] = (uzelx[N/2]+uzelx[N/2+1])/2.0;
//      System.out.println("uzelxN/2 = "+ uzelx[N/2]);

        for (int i=1;i<N+1;i++){
            h[i] = uzelx[i]-uzelx[i-1];
        }

        return h;

    }

    public static double[] setkaBakhvalova(double epsilon, int N){
//        System.out.println("Bakhvalov mesh");
        double uzelx[] = new double[N+1];
        double sigma1;
        double[] h = new double[N+1];
        sigma1 = Math.min(0.5,(-4.)*epsilon*Math.log(epsilon)/1.);
//        sigma2 = min(0.5,(-2.)*epsilon*log(epsilon)/2.);
        if (epsilon>Math.exp(-1)) {
            sigma1 = 0.5;
//            sigma2 = 0.5;
        }


        if (sigma1==0.5){
            uzelx[0]=0.;
            for(int i=1;i<=N;i++){
                h[i]=1./N;
                uzelx[i]=uzelx[i-1]+h[i];
            }
        }


        if (sigma1<0.5)
        {
            uzelx[0] = 0.;
            for (int i=1;i<=N/2;i++)
            {
                uzelx[i] = (-4.)/1.0*epsilon*Math.log(1.-2.*(1.-epsilon)*i/N);
                //cout<<uzelx[i]<<endl;
            }
            for (int i=N/2;i<=N;i++)
            {
                uzelx[i] = sigma1 + (2.*i/N - 1.)*(1.-sigma1);
            }


//            uzelx[N/2-1] = (uzelx[N/2-1] + uzelx[N/2])/2.0;
//            uzelx[N/2] = (uzelx[N/2]+uzelx[N/2+1])/2.0;

//            uzelx[N/2-1] = (uzelx[N/2-1] + uzelx[N/2])/2.0;
//            uzelx[N/2] = (uzelx[N/2]+uzelx[N/2+1])/2.0;


            for(int i=1;i<=N;i++)
            {
                h[i]=uzelx[i]-uzelx[i-1];
            }
        }


        return h;
    }

    public static void findPoints(int N, double[] h, double[] uzel){
        uzel[0]= 0.;
        for (int i=1;i<N+1;i++)
        {
            uzel[i] = uzel[i-1]+h[i];
//            System.out.println(uzel[i]);
        }
    }


    public static void findFunction(int N, double[] f, double[] uzel, Function<Double, Double> function){
        for (int i=0;i<N+1;i++)
        {
            f[i] = function.apply(uzel[i]);
        }
    }

    public static void findCoeff(double[] h, double[] f, double[] A, double[] B, double[] C, double[] F, int N){
        for (int i = 1; i < N; i++) {
            A[i] = h[i];
            C[i] = (-2.)*(h[i] + h[i+1]);
            B[i] = h[i+1];
            F[i] = (-6.) * ( (f[i + 1] - f[i])/h[i+1] - (f[i] - f[i - 1])/h[i] );

        }
    }

    public static void findS(int L, int N, double[] uzel, double[] x, double[] f, double[] S, double[] h, double[] hh, double[] c){

        x[0]=0.;

        for(int k=1;k<L+1;k++){
            x[k]=x[k-1]+hh[k];
        }

        for(int i=0;i<N+1;i++){
            for (int k=0;k<L;k++){
                // if ( ((x[k]>=uzel[i])&&(x[k]<uzel[i+1]))||(x[k]==1.) )
                if  ((x[k]>=uzel[i])&&(x[k]<uzel[i+1]))
                {

//                    S[k] =
                    S[k] = (uzel[i+1] - x[k])*(uzel[i+1] - x[k])*(uzel[i+1] - x[k])*c[i]/(6.*h[i+1]) + (x[k] - uzel[i])*(x[k] - uzel[i])*(x[k] - uzel[i])*c[i+1]/(6.*h[i+1]) + (f[i+1]/h[i+1] - c[i+1]*h[i+1]/6.)*(x[k] - uzel[i]) + (f[i]/h[i+1] - c[i]*h[i+1]/6.)*(uzel[i+1] - x[k]);

                }
            }
        }

        S[L] = (uzel[N] - x[L])*(uzel[N] - x[L])*(uzel[N] - x[L])*c[N]/(6.*h[N]) + (x[L] - uzel[N-1])*(x[L] - uzel[N-1])*(x[L] - uzel[N-1])*c[N]/(6.*h[N]) + (f[N]/h[N] - c[N]*h[N]/6.)*(x[L] - uzel[N-1]) + (f[N-1]/h[N] - c[N-1]*h[N]/6.)*(uzel[N] - x[L]);

    }

    public static void pogreshnost(double[] res, int L, double[] x, double[] S, int N, Function<Double, Double> function, List<Double> norm, double[] uzel, double[] f){
        AkimaSplineInterpolator interpolator = new AkimaSplineInterpolator();
//        NaturalSplineInterpolator interpolator = new NaturalSplineInterpolator();
        PolynomialSplineFunction spline = interpolator.interpolate(uzel, f);
        double norma = 0.;
        for(int i=0;i<L+1;i++)
        {
            res[i] = Math.abs(function.apply(x[i]) - S[i]);
//            res[i] = Math.abs(function.apply(x[i]) - spline.value(x[i]));
//            System.out.println("res["+i+"] = "+res[i]);
        }

        //норма погрешности:
        for (int i=0;i<L+1;i++)
        {
            if(res[i]>norma) norma=res[i];
//            System.out.println(i+" = "+norma);
        }
        norm.add(norma);
    }



    public static double[] cubicalSpline(int N, Function<Double, Double> function, double epsilon, List<Double> norm){
        final int L = (5*N);
        double[] x = new double[L+1];

        double[]A = new double[N+1];
        double[]B = new double[N+1];
        double[]C = new double[N+1];
        double[]F = new double[N+1];
        double[]f = new double[N+1];
        double[]uzel = new double[N+1];
        double[]c = new double[N+1];

        double[]res = new double[L+1];
        double[]S = new double[L+1];
        double[] hRavnom = ravnomSetka(N);
        double[] hRavnomGust = ravnomSetka(L);
//        double h = 1./N;
//        double epsilon = 10.e-2;
        double[] hSetkaShishkina = setkaShishkina(epsilon, N);
        double[] hSetkaShishkinaGust = setkaShishkina(epsilon, L);
        double[] hSetkaBakhvalova = setkaBakhvalova(epsilon, N);
        double[] hSetkaBakhvalovaGust = setkaBakhvalova(epsilon, L);

        int iteration = 0;
//        double[]h = hSetkaBakhvalova;
//        double[]hGust = hSetkaBakhvalovaGust;
//        double[]h = hSetkaShishkina;
//        double[]hGust = hSetkaShishkinaGust;
        double[]h = hRavnom;
        double[]hGust = hRavnomGust;

        findPoints(N, h, uzel);
        findFunction(N, f, uzel, function);
        findCoeff(h, f, A, B, C, F, N);
        progonka(A, B, C, F, N, c, f, h, epsilon, iteration);
        findS(L, N, uzel, x, f, S, h, hGust, c);
        pogreshnost(res, L, x, S, N, function, norm, uzel, f);
//        System.out.println("iteration = "+iteration);
        return S;
    }




    public static void main(String[] args){
        Map<Integer, List> map = new HashMap<>();
        double epsilon = 1.;
        System.out.println("Epsilon = "+epsilon);
        double a = 0.;
        for (int i=16; i<=512; i=i*2){
            List<Double> norm = new ArrayList<>();
//            cubicalSpline(i, x->x*x*x*x*x, epsilon, norm);
//            cubicalSpline(i, x->Math.cos(x+0.5), epsilon, norm);
//            cubicalSpline(i, x->Math.cos(Math.PI*x/4.), epsilon, norm);
//            cubicalSpline(i, x->Math.sin(Math.PI*x/2.) + Math.exp(-x/epsilon), epsilon, norm);
//            cubicalSpline(i, x->Math.sin(x)*Math.sin(x), epsilon, norm);
//            cubicalSpline(i, x->Math.sin(x), epsilon, norm);
//            cubicalSpline(i, x->Math.sin(x)*Math.sin(x)*Math.sin(x), epsilon, norm);//4
//            cubicalSpline(i, x->Math.cos(Math.PI*x/2.) + Math.exp(-x/epsilon), epsilon, norm);
            cubicalSpline(i, x->Math.sin(Math.PI*x/2.) + Math.exp(-x/epsilon), epsilon, norm);
//            cubicalSpline(i, x->Math.cos(Math.PI*x/1.) + Math.exp(-x/epsilon), epsilon, norm);
            double four = a/norm.get(0);
            a = norm.get(0);
            System.out.println("Порядок точности log2 (||"+i/2.+"||/||"+i+"||) = "+Math.log10(four)/Math.log10(2.));
            System.out.println("_______________________________________________________________");
            System.out.println("||"+i+"|| = "+ norm);

        }
    }
}
