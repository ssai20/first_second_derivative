import org.apache.commons.math3.analysis.interpolation.AkimaSplineInterpolator;
import org.apache.commons.math3.analysis.interpolation.BicubicSplineInterpolatingFunction;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

import java.util.*;
import java.util.function.BiFunction;
import java.util.function.Function;

public class BiCubicSpline {

        public static void progonka(double[] A, double[] B, double[] C, double[] F, int N, double[] c, double[] f, double[] h, double epsilon, Function<Double, Double> secondDerivation0, Function<Double, Double> secondDerivationN, Function<Double, Double> betta1, double xory){
            double[]alpha = new double[N+1];
            double[]betta = new double[N+1];

//            betta[1] = 3.0/h[1]*((f[0]-f[1])/h[1]-betta1.apply(epsilon));
//            alpha[1] = -0.5;


//            betta[1] = -Math.PI*Math.PI/4.+1.0/(epsilon*epsilon);//2
//            betta[1] = secondDerivation0.apply(epsilon);
//            System.out.println(xory);
            alpha[1] = 0.0;
            betta[1] = betta1.apply(xory);
            for(int i=1;i<N;i++)
            {
                alpha[i+1]=B[i]/(C[i]-alpha[i]*A[i]);

                betta[i+1]=(A[i]*betta[i]+F[i])/(C[i]-alpha[i]*A[i]);
            }

//            c[0]=-Math.PI*Math.PI/4.0+1.0/(epsilon*epsilon); //вторая производная cos(Pi*x/4)+epsilon в 0
//            c[N]=Math.exp(-1./epsilon)/(epsilon*epsilon);//вторая производная cos(Pi*x/4)+epsilon в 1

//            c[0] = secondDerivation0.apply(xory);
//            c[N] = secondDerivationN.apply(xory);

//            c[0] = 0.0;
//            c[N] = 0.0;

//        c[0] = 0.;//sin^3
//        c[N] = 3.*(2.*Math.sin(1.)*Math.cos(1.)*Math.cos(1.) - Math.sin(1.)*Math.sin(1.)*Math.sin(1.));//sin^3
//        System.out.println(3.*(2.*Math.sin(1.)*Math.cos(1.)*Math.cos(1.) - Math.sin(1.)*Math.sin(1.)*Math.sin(1.)));


            double proizv1 = secondDerivation0.apply(epsilon);
            double proizv2 = secondDerivationN.apply(epsilon);


//        double proizv1 = -1.0/epsilon; //производная Math.cos(Math.PI*x/2.) + Math.exp(-x/epsilon) в 0
//        double proizv2 = -Math.PI/2.- Math.exp(-1.0/epsilon)/epsilon; //производная Math.cos(Math.PI*x/2.) + Math.exp(-x/epsilon) в 1
//
        c[0] = 6.*alpha[1]/(h[1]*(2.*alpha[1]+1.)) * ((h[1]*betta[1])/(6.*alpha[1]) + (f[1]-f[0])/h[1] - proizv1);
        c[N] = 6.*((f[N-1]-f[N])/h[N] - h[N]*betta[N]/6. + proizv2)/(h[N]*(2.+alpha[N]));


            for (int i=N-1;i>0;i--)
            {
                c[i]=alpha[i+1]*c[i+1] + betta[i+1];
//                System.out.println(c[i]);
            }
        }

        public static  double[] ravnomSetka(int N){
            double[] h = new double[N+1];
            for (int i=1;i<N+1;i++){
                h[i] = 1.0/N;
            }
            return h;
        }


        public static List<double[]> setkaShishkina(double epsilon, int N, int M){
//      System.out.println("Shishkin mesh");
//            double sigma1;
            double[] h = new double[N+1];
            double[] t = new double[M+1];
            double sigma1 = Math.min(0.5, 4.*epsilon*Math.log(N*1.0)/1.);
            double sigma2 = Math.min(0.5, 4.*epsilon*Math.log(M*1.0)/2.);

            for(int i=1;i<=N;i++){
                if ((i<=N/2)&&(i>=1)) {
                    h[i] = 2.*sigma1/N;
                    t[i] = 2.*sigma2/N;
                }
                if ((i<=N)&&(i>N/2)) {
                    h[i] = 2.*(1.-sigma1)/N;
                    t[i] = 2.*(1.-sigma2)/N;
                }

//                System.out.println(h[i]);
            }

//            double uzelx[] = new double[N+1];
//            uzelx[0]=0.;
////            uzely[0]=0.;
//            for (int i=1;i<N+1;i++)
//            {
//                uzelx[i]=uzelx[i-1]+h[i];
//
//            }
            List<double[]> list = new ArrayList<>();
            list.add(h);
            list.add(t);

            return list;

        }

    public static List<double[]> setkaShishkinaModern(double epsilon, int N, int M){
//      System.out.println("Shishkin mesh");
//            double sigma1;
        double[] h = new double[N+1];
        double[] t = new double[M+1];
        double sigma1 = Math.min(0.5, 4.*epsilon*Math.log(N*1.0)/1.);
        double sigma2 = Math.min(0.5, 4.*epsilon*Math.log(M*1.0)/2.);

        for(int i=1;i<=N;i++){
            if ((i<=N/2)&&(i>=1)) {
                h[i] = 2.*sigma1/N;
                t[i] = 2.*sigma2/N;
            }
            if ((i<=N)&&(i>N/2)) {
                h[i] = 2.*(1.-sigma1)/N;
                t[i] = 2.*(1.-sigma2)/N;
            }

        }

            double uzelx[] = new double[N+1];
            uzelx[0]=0.0;
            for (int i=1;i<N+1;i++)
            {
                uzelx[i]=uzelx[i-1]+h[i];
            }
            uzelx[N/2] = (uzelx[N/2]+uzelx[N/2+1])/2.0;      // модернизация со сдвигом вправо
            for (int i=1;i<N+1;i++){
                h[i] = uzelx[i]-uzelx[i-1];
            }

        double uzely[] = new double[N+1];
        uzely[0]=0.0;
        for (int i=1;i<N+1;i++)
        {
            uzely[i]=uzely[i-1]+t[i];
        }
        uzely[N/2] = (uzely[N/2]+uzely[N/2+1])/2.0;
        for (int i=1;i<N+1;i++){
            t[i] = uzely[i]-uzely[i-1];
        }




        List<double[]> list = new ArrayList<>();
        list.add(h);
        list.add(t);

        return list;

    }


    public static List<double[]> setkaBakhvalova(double epsilon, int N, int M){
//        System.out.println("Bakhvalov mesh");
            double uzelx[] = new double[N+1];
            double uzely[] = new double[M+1];
//            double sigma1;
            double[] h = new double[N+1];
            double[] t = new double[M+1];
            double sigma1 = Math.min(0.5,(-4.)*epsilon*Math.log(epsilon)/1.);
            double sigma2 = Math.min(0.5,(-4.)*epsilon*Math.log(epsilon)/2.);
            if (epsilon>Math.exp(-1)) {
                sigma1 = 0.5;
                sigma2 = 0.5;
            }


            if (sigma1==0.5){
                uzelx[0]=0.;
                for(int i=1;i<=N;i++){
                    h[i]=1./N;
                    uzelx[i]=uzelx[i-1]+h[i];
                }
            }

            if (sigma2==0.5){
                uzely[0]=0.;
                for(int i=1;i<=M;i++){
                    t[i]=1./M;
                    uzely[i]=uzely[i-1]+t[i];
                }
            }


            if (sigma1<0.5)
            {
                uzelx[0] = 0.;
                for (int i=1;i<=N/2;i++)
                {
                    uzelx[i] = (-4.0)/1.0*epsilon*Math.log(1.-2.*(1.-epsilon)*i/N);
                    //cout<<uzelx[i]<<endl;
                }
                for (int i=N/2;i<=N;i++)
                {
                    uzelx[i] = sigma1 + (2.*i/N - 1.)*(1.-sigma1);
                }
                for(int i=1;i<=N;i++)
                {
                    h[i]=uzelx[i]-uzelx[i-1];
                }
            }

            if (sigma2<0.5)
            {
                uzely[0] = 0.;
                for (int i=1;i<=M/2;i++)
                {
                    uzely[i] = (-4.0)/2.0*epsilon*Math.log(1.-2.*(1.-epsilon)*i/M);
                    //cout<<uzelx[i]<<endl;
                }
                for (int i=M/2;i<=M;i++)
                {
                    uzely[i] = sigma2 + (2.*i/M - 1.)*(1.-sigma2);
                }
                for(int i=1;i<=M;i++)
                {
                    t[i]=uzely[i]-uzely[i-1];
                }
            }

            List<double[]> list = new ArrayList<>();
            list.add(h);
            list.add(t);
            return list;
        }


    public static List<double[]> setkaBakhvalovaModern(double epsilon, int N, int M){
//        System.out.println("Bakhvalov mesh");
        double uzelx[] = new double[N+1];
        double uzely[] = new double[M+1];
//            double sigma1;
        double[] h = new double[N+1];
        double[] t = new double[M+1];
        double sigma1 = Math.min(0.5,(-4.)*epsilon*Math.log(epsilon)/1.);
        double sigma2 = Math.min(0.5,(-4.)*epsilon*Math.log(epsilon)/2.);
        if (epsilon>Math.exp(-1)) {
            sigma1 = 0.5;
            sigma2 = 0.5;
        }


        if (sigma1==0.5){
            uzelx[0]=0.;
            for(int i=1;i<=N;i++){
                h[i]=1./N;
                uzelx[i]=uzelx[i-1]+h[i];
            }
            uzelx[N/2-1] = (uzelx[N/2-1] + uzelx[N/2])/2.0;
            uzelx[N/2] = (uzelx[N/2]+uzelx[N/2-1])/2.0;
//            \bar{x_i}= (x_i+x_{i-1})/2   i= N/2-1, i=N/2
        }

        if (sigma2==0.5){
            uzely[0]=0.;
            for(int i=1;i<=M;i++){
                t[i]=1./M;
                uzely[i]=uzely[i-1]+t[i];
            }
            uzely[M/2-1] = (uzely[M/2-1] + uzely[M/2])/2.0;
            uzely[M/2] = (uzely[M/2]+uzely[M/2+1])/2.0;
        }


        if (sigma1<0.5)
        {
            uzelx[0] = 0.;
            for (int i=1;i<=N/2;i++)
            {
                uzelx[i] = (-4.0)/1.0*epsilon*Math.log(1.-2.*(1.-epsilon)*i/N);
                //cout<<uzelx[i]<<endl;
            }
            for (int i=N/2;i<=N;i++)
            {
                uzelx[i] = sigma1 + (2.*i/N - 1.)*(1.-sigma1);
            }
            uzelx[N/2-1] = (uzelx[N/2-1] + uzelx[N/2])/2.0;
            uzelx[N/2] = (uzelx[N/2]+uzelx[N/2+1])/2.0;
            for(int i=1;i<=N;i++)
            {
                h[i]=uzelx[i]-uzelx[i-1];
            }
        }

        if (sigma2<0.5)
        {
            uzely[0] = 0.;
            for (int i=1;i<=M/2;i++)
            {
                uzely[i] = (-4.0)/2.0*epsilon*Math.log(1.-2.*(1.-epsilon)*i/M);
                //cout<<uzelx[i]<<endl;
            }
            for (int i=M/2;i<=M;i++)
            {
                uzely[i] = sigma2 + (2.*i/M - 1.)*(1.-sigma2);
            }
            uzely[M/2-1] = (uzely[M/2-1] + uzely[M/2])/2.0;
            uzely[M/2] = (uzely[M/2]+uzely[M/2+1])/2.0;
            for(int i=1;i<=M;i++)
            {
                t[i]=uzely[i]-uzely[i-1];
            }
        }



//        double uzelx[] = new double[N+1];
//        uzelx[0]=0.0;
//        for (int i=1;i<N+1;i++)
//        {
//            uzelx[i]=uzelx[i-1]+h[i];
//        }
//        uzelx[N/2-1] = (uzelx[N/2-1] + uzelx[N/2])/2.0;
//        uzelx[N/2] = (uzelx[N/2]+uzelx[N/2+1])/2.0;
//        for (int i=1;i<N+1;i++){
//            h[i] = uzelx[i]-uzelx[i-1];
//        }
//
//        double uzely[] = new double[N+1];
//        uzely[0]=0.0;
//        for (int i=1;i<N+1;i++)
//        {
//            uzely[i]=uzely[i-1]+t[i];
//        }
//        uzely[N/2-1] = (uzely[N/2-1] + uzely[N/2])/2.0;
//        uzely[N/2] = (uzely[N/2]+uzely[N/2+1])/2.0;
//        for (int i=1;i<N+1;i++){
//            t[i] = uzely[i]-uzely[i-1];
//        }



        List<double[]> list = new ArrayList<>();
        list.add(h);
        list.add(t);
        return list;
    }



    public static void findPoints(int N, double[] h, double[] uzel){
            uzel[0] = 0.;
            for (int i=1;i<N+1;i++)
            {
                uzel[i] = uzel[i-1]+h[i];
//            System.out.println(uzel[i]);
            }
        }


    public static void findPointsShishkinaModern(int N, double[] h, double[] uzel){
        uzel[0] = 0.;
        for (int i=1;i<N+1;i++)
        {
            uzel[i] = uzel[i-1]+h[i];
        }
        uzel[N/2] = (uzel[N/2]+uzel[N/2+1])/2.0;
    }

    public static void findPointsBakhvalovaModern(int N, double[] h, double[] uzel){
        uzel[0] = 0.;
        for (int i=1;i<N+1;i++)
        {
            uzel[i] = uzel[i-1]+h[i];
        }
        uzel[N/2-1] = (uzel[N/2-1] + uzel[N/2])/2.0;
        uzel[N/2] = (uzel[N/2]+uzel[N/2+1])/2.0;
    }

        public static void findFunction(int N, double[] f, double[] uzel, Function<Double, Double> function){
            for (int i=0;i<N+1;i++)
            {
                f[i] = function.apply(uzel[i]);
            }
        }

    public static void findFunction2(int N, int M, double[][] f2, double[] uzelX, double[] uzelY, BiFunction<Double, Double, Double> function2){
        for (int i=0;i<N+1;i++)
        {
            for (int j=0;j<M+1;j++)
            {
                f2[i][j] = function2.apply(uzelX[i],uzelY[j]);
//                System.out.println(f2[i][j]);
            }
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

            for(int i=0;i<N+1;i++){//проверить
                for (int k=0;k<L;k++){ //проверить этот момент
                    if  ((x[k]>=uzel[i])&&(x[k]<uzel[i+1]))
                    {
                        S[k] = (uzel[i+1] - x[k])*(uzel[i+1] - x[k])*(uzel[i+1] - x[k])*c[i]/(6.*h[i+1]) + (x[k] - uzel[i])*(x[k] - uzel[i])*(x[k] - uzel[i])*c[i+1]/(6.*h[i+1]) + (f[i+1]/h[i+1] - c[i+1]*h[i+1]/6.)*(x[k] - uzel[i]) + (f[i]/h[i+1] - c[i]*h[i+1]/6.)*(uzel[i+1] - x[k]);
                    }
                }
            }

            S[L] = (uzel[N] - x[L])*(uzel[N] - x[L])*(uzel[N] - x[L])*c[N]/(6.*h[N]) + (x[L] - uzel[N-1])*(x[L] - uzel[N-1])*(x[L] - uzel[N-1])*c[N]/(6.*h[N]) + (f[N]/h[N] - c[N]*h[N]/6.)*(x[L] - uzel[N-1]) + (f[N-1]/h[N] - c[N-1]*h[N]/6.)*(uzel[N] - x[L]);

        }


        public static void pogreshnost(double[] res, int L, double[] x, double[] S, int N, Function<Double, Double> function, List<Double> norm, double[] uzel, double[] f){
//            AkimaSplineInterpolator interpolator = new AkimaSplineInterpolator();
////        NaturalSplineInterpolator interpolator = new NaturalSplineInterpolator();
//            PolynomialSplineFunction spline = interpolator.interpolate(uzel, f);
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

    public static void pogreshnost2(double[][] res2, int L, int K, double[] x, double[] y, double[][] U, BiFunction<Double, Double, Double> function2, List<Double> norm){

//        BicubicSplineInterpolatingFunction interpolator = new BicubicSplineInterpolatingFunction(uzelX, uzelY, function2);
        double norma = 0.0;
        for(int i=0;i<L+1;i++)
        {
//            System.out.println("x["+i+"] = "+x[i]);
            for(int j=0;j<K+1;j++) {
                res2[i][j] = Math.abs(function2.apply(x[i], y[j]) - U[i][j]);

//                System.out.print(function2.apply(x[i], y[j]));
//                System.out.print("  ["+i+"]"+"["+j+"]   ");
//                System.out.print(U[i][j]);
//                System.out.println("  "+" = "+res2[i][j]);


//            res2[i] = Math.abs(function2.apply(x[i],y[j]) - interpolator.value(x[i], y[j]));

//                System.out.println("y["+j+"] = "+y[j]);
//            System.out.println("res["+i+"] = "+res[i]);
            }
        }

        //норма погрешности:
        double ii = 0,jj=0;
        for(int i=0;i<L+1;i++)
        {
            for(int j=0;j<K+1;j++) {
                if (res2[i][j] > norma) {
                    norma = res2[i][j];
                    ii=i;
                    jj=j;
                }
//            System.out.println(i+" = "+norma);
            }
        }
        norm.add(norma);
        norm.add(ii);
        norm.add(jj);
    }



        public static double[] cubicalSpline(int N, Function<Double, Double> function, double epsilon, List<Double> norm, Function<Double, Double> secondDerivation0, Function<Double, Double> secondDerivationN, String setka){
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

            double[] hSetkaShishkina = setkaShishkina(epsilon, N, N).get(0);
            double[] hSetkaShishkinaGust = setkaShishkina(epsilon, L, L).get(0);
            double[] hSetkaBakhvalova = setkaBakhvalova(epsilon, N, N).get(0);
            double[] hSetkaBakhvalovaGust = setkaBakhvalova(epsilon, L, L).get(0);

            double[] h = new double[N+1];
            double[]hGust = new double[L+1];

            if (setka == "b") {
                h = hSetkaBakhvalova;
                hGust = hSetkaBakhvalovaGust;
            }
            if (setka == "s") {
                h = hSetkaShishkina;
                hGust = hSetkaShishkinaGust;
            }
            if (setka == "r") {
                h = hRavnom;
                hGust = hRavnomGust;
            }

            findPoints(N, h, uzel);
            findFunction(N, f, uzel, function);
            findCoeff(h, f, A, B, C, F, N);
            progonka(A, B, C, F, N, c, f, h, epsilon, secondDerivation0, secondDerivationN, secondDerivation0, 0);

            findS(L, N, uzel, x, f, S, h, hGust, c);
            pogreshnost(res, L, x, S, N, function, norm, uzel, f);
            return S;
        }

        public static void findGl(double[] G1i, double[] G2i, double[] G3i, double uzelXi, Function<Double, Double> G1i0, Function<Double, Double> G3i0, int N, int M, double[] f2i, double[] t, double[] Moment, double epsilon){
//            double[][] G0 = new double[N+1][N+1];
//            double[][] G1 = new double[N+1][N+1];
//            double[][] G2 = new double[N+1][N+1];
//            double[][] G3 = new double[N+1][N+1];
//            G0 = f2;
            for(int j=1;j<M+1;j++)
                {
                    G1i[j] = t[j]*Moment[j]/2.0 + f2i[j]/t[j] - Moment[j]*t[j]/6.0-f2i[j-1]/t[j] + Moment[j-1]*t[j]/6.0;
//                    System.out.println(G1i[j]);
                    G2i[j] = Moment[j];
//                    System.out.println(G2i[j]);
                    G3i[j] = (Moment[j] - Moment[j-1])/t[j];
                }
//            System.out.println(uzelXi);
            G1i[0] = G1i0.apply(uzelXi);
            G2i[0] = Moment[0];
            G3i[0] = G3i0.apply(uzelXi);




//            G1i[0] = -1.0/epsilon;
//            G2i[0] = Moment[0];
//            G3i[0] = -1.0/(epsilon*epsilon*epsilon);

        }
    public static double[][] bicubicalSpline(int N, int M, int L, int K, double epsilon, BiFunction<Double, Double, Double> biFunction, Function<Double, Double> derivationY0, Function<Double, Double> derivationYN, Function<Double, Double> betta1Y, Function<Double, Double> nullDerivationX0, Function<Double, Double> nullDerivationXN, Function<Double, Double> nullBetta1X, Function<Double, Double> firstDerivationX0, Function<Double, Double> firstDerivationXN, Function<Double, Double> firstBetta1X, Function<Double, Double> secondDerivationX0, Function<Double, Double> secondDerivationXN, Function<Double, Double> secondBetta1X, Function<Double, Double> thirdDerivationX0, Function<Double, Double> thirdDerivationXN, Function<Double, Double> thirdBetta1X, List<Double> norm, String setkaX, String setkaY, Function<Double, Double> G1i0, Function<Double, Double> G3i0, Function<Double, Double> a0y10, Function<Double, Double> a0y30) {
//        final int K = 3 * M;
//        final int L = 3 * N;
//        final double epsilon = 1.e-2;

        double[][] u = new double[L + 1][K + 1];

        double[][] a00 = new double[N + 1][M + 1];
        double[][] a01 = new double[N + 1][M + 1];
        double[][] a02 = new double[N + 1][M + 1];
        double[][] a03 = new double[N + 1][M + 1];
        double[][] a10 = new double[N + 1][M + 1];
        double[][] a11 = new double[N + 1][M + 1];
        double[][] a12 = new double[N + 1][M + 1];
        double[][] a13 = new double[N + 1][M + 1];
        double[][] a20 = new double[N + 1][M + 1];
        double[][] a21 = new double[N + 1][M + 1];
        double[][] a22 = new double[N + 1][M + 1];
        double[][] a23 = new double[N + 1][M + 1];
        double[][] a30 = new double[N + 1][M + 1];
        double[][] a31 = new double[N + 1][M + 1];
        double[][] a32 = new double[N + 1][M + 1];
        double[][] a33 = new double[N + 1][M + 1];

        double[] uzelX = new double[N + 1];
        double[] uzelY = new double[M + 1];
        double[] t = new double[M + 1];
        double[] h = new double[N + 1];
        double[] tGust = new double[K + 1];
        double[] hGust = new double[L + 1];

        double[] A = new double[N + 1];
        double[] B = new double[N + 1];
        double[] C = new double[N + 1];
        double[] F = new double[N + 1];

        double[][] f2 = new double[N + 1][M + 1];

        double[] Moment = new double[M + 1];
        double[] M0 = new double[N + 1];
        double[] M1 = new double[N + 1];
        double[] M2 = new double[N + 1];
        double[] M3 = new double[N + 1];

        double[][] G0 = new double[N + 1][M + 1];
        double[][] G1 = new double[N + 1][M + 1];
        double[][] G2 = new double[N + 1][M + 1];
        double[][] G3 = new double[N + 1][M + 1];

        double[][] res2 = new double[L + 1][K + 1];

        double[] x = new double[L+1];
        double[] y = new double[K+1];

        double[] hRavnom = ravnomSetka(N);
        double[] hRavnomGust = ravnomSetka(L);
        double[] hSetkaShishkina = setkaShishkina(epsilon, N, M).get(0);
        double[] hSetkaShishkinaGust = setkaShishkina(epsilon, L, K).get(0);
        double[] hSetkaBakhvalova = setkaBakhvalova(epsilon, N, M).get(0);
        double[] hSetkaBakhvalovaGust = setkaBakhvalova(epsilon, L, K).get(0);

        double[] hSetkaShishkinaModern = setkaShishkinaModern(epsilon, N, M).get(0);
        double[] hSetkaShishkinaGustModern = setkaShishkinaModern(epsilon, L, K).get(0);
        double[] hSetkaBakhvalovaModern = setkaBakhvalovaModern(epsilon, N, M).get(0);
        double[] hSetkaBakhvalovaGustModern = setkaBakhvalovaModern(epsilon, L, K).get(0);

        double[] tRavnom = ravnomSetka(M);
        double[] tRavnomGust = ravnomSetka(K);
        double[] tSetkaShishkina = setkaShishkina(epsilon, N, M).get(1);
        double[] tSetkaShishkinaGust = setkaShishkina(epsilon, L, K).get(1);
        double[] tSetkaBakhvalova = setkaBakhvalova(epsilon, N, M).get(1);
        double[] tSetkaBakhvalovaGust = setkaBakhvalova(epsilon, L, K).get(1);

        double[] tSetkaShishkinaModern = setkaShishkinaModern(epsilon, N, M).get(1);
        double[] tSetkaShishkinaGustModern = setkaShishkinaModern(epsilon, L, K).get(1);
        double[] tSetkaBakhvalovaModern = setkaBakhvalovaModern(epsilon, N, M).get(1);
        double[] tSetkaBakhvalovaGustModern = setkaBakhvalovaModern(epsilon, L, K).get(1);

        if (setkaX == "b") {
            h = hSetkaBakhvalova;
            hGust = hSetkaBakhvalovaGust;
        }

        if (setkaX == "bm") {
            h = hSetkaBakhvalovaModern;
            hGust = hSetkaBakhvalovaGustModern;
        }

        if (setkaX == "s") {
            h = hSetkaShishkina;
            hGust =hSetkaShishkinaGust;
        }

        if (setkaX == "sm") {
            h = hSetkaShishkinaModern;
            hGust =hSetkaShishkinaGustModern;
        }

        if (setkaX == "r") {
            h = hRavnom;
            hGust = hRavnomGust;
        }
        if (setkaY == "b") {
            t = tSetkaBakhvalova;
            tGust = tSetkaBakhvalovaGust;
        }

        if (setkaY == "bm") {
            t = tSetkaBakhvalovaModern;
            tGust = tSetkaBakhvalovaGustModern;
        }



        if (setkaY == "s") {
            t = tSetkaShishkina;
            tGust = tSetkaShishkinaGust;
        }

        if (setkaY == "sm") {
            t = tSetkaShishkinaModern;
            tGust = tSetkaShishkinaGustModern;
        }

        if (setkaY == "r") {
            t = tRavnom;
            tGust = tRavnomGust;
        }

        findPoints(N, h, uzelX);
        findPoints(M, t, uzelY);
        findPoints(L, hGust, x);
        findPoints(K, tGust, y);

//        findPointsShishkinaModern(N, h, uzelX);
//        findPointsShishkinaModern(M, t, uzelY);
//        findPointsShishkinaModern(L, hGust, x);
//        findPointsShishkinaModern(K, tGust, y);

//        findPointsBakhvalovaModern(N, h, uzelX);
//        findPointsBakhvalovaModern(M, t, uzelY);
//        findPointsBakhvalovaModern(L, hGust, x);
//        findPointsBakhvalovaModern(K, tGust, y);


//        for(int i=0;i<L+1;i++)
//        {
//            System.out.println("x["+i+"] = "+x[i]);
//        }
//        for(int j=0;j<K+1;j++)
//        {
//            System.out.println("y["+j+"] = "+x[j]);
//        }

        findFunction2(N, M, f2, uzelX, uzelY, biFunction);

//        G0 = f2;
        for (int i = 0; i < N + 1; i++) {
            findCoeff(t, f2[i], A, B, C, F, M);
            progonka(A, B, C, F, M, Moment, f2[i], t, epsilon, derivationY0, derivationYN, betta1Y, uzelX[i]);
//            progonka(A, B, C, F, M, Moment, f2[i], t, epsilon, eps -> -1.0 / (eps), eps -> -Math.exp(-1.0 / eps) / (eps), eps -> 1.0 / (eps * eps));
            findGl(G1[i], G2[i], G3[i], uzelX[i], G1i0, G3i0, N, M, f2[i], t, Moment, epsilon);
        }
//        for (int i=0;i<N+1;i++)
//        {
//            for(int j=0;j<M+1;j++)
//            {
//                System.out.println(G3[i][j]);
//            }
//        }

        G0 = f2;
        for (int j = 0; j < M + 1; j++) {
            double[] column = new double[N+1];
            double[] column1 = new double[N+1];
            for (int i=0;i<N+1;i++)
            {
                column[i] = G0[i][j];
            }
            findCoeff(h, column, A, B, C, F, N);
            progonka(A, B, C, F, N, M0, column, h, epsilon, nullDerivationX0, nullDerivationXN, nullBetta1X, uzelY[j]);
//            progonka(A, B, C, F, N, M0, column, h, epsilon, eps ->0.0, eps -> -Math.PI / 2.0, eps -> -Math.PI * Math.PI / 4.0);
            for (int i=0;i<N+1;i++)
            {
                column[i] = G1[i][j];
//                System.out.println("G1 = "+G1[i][j]);
            }
            findCoeff(h, column, A, B, C, F, N);
            progonka(A, B, C, F, N, M1, column, h, epsilon, firstDerivationX0, firstDerivationXN, firstBetta1X, uzelY[j]);
//            progonka(A, B, C, F, N, M1, column, h, epsilon, eps->0.0, eps->0.0, eps->0.0);
//            progonka(A, B, C, F, N, M0, G0[j], h, epsilon, eps ->0.0, eps -> -Math.PI / 2.0, eps -> -Math.PI * Math.PI / 4.0);
            for (int i=0;i<N+1;i++)
            {
                column[i] = G2[i][j];
//                System.out.println(G2[i][j]);
            }
            findCoeff(h, column, A, B, C, F, N);
            progonka(A, B, C, F, N, M2, column, h, epsilon, secondDerivationX0, secondDerivationXN, secondBetta1X, uzelY[j]);
//            progonka(A, B, C, F, N, M1, column, h, epsilon, eps->0.0, eps->0.0, eps->0.0);
//            progonka(A, B, C, F, N, M0, G0[j], h, epsilon, eps ->0.0, eps -> -Math.PI / 2.0, eps -> -Math.PI * Math.PI / 4.0);
            for (int i=0;i<N+1;i++)
            {
                column[i] = G3[i][j];
            }
            findCoeff(h, column, A, B, C, F, N);
//            progonka(A, B, C, F, N, M3, column, h, epsilon, thirdDerivationX0, thirdDerivationXN, thirdBetta1X, uzelY[j]);
//            progonka(A, B, C, F, N, M1, column, h, epsilon, eps->0.0, eps->0.0, eps->0.0);

//            progonka(A, B, C, F, N, M0, G0[j], h, epsilon, eps ->0.0, eps -> -Math.PI / 2.0, eps -> -Math.PI * Math.PI / 4.0);
//            System.out.println("M["+j+"] = "+M3[0][j]);
        }

        for (int i = 1; i < N+1; i++) {
            for (int j = 0; j < M+1; j++) {
                a00[i][j] = G0[i][j];
                a01[i][j] = G1[i][j];
                a02[i][j] = G2[i][j];
                a03[i][j] = G3[i][j];

                a10[i][j] = h[i] * M0[i] / 2.0 + G0[i][j] / h[i] - M0[i] * h[i] / 6.0 - G0[i - 1][j] / h[i] + M0[i - 1] * h[i] / 6.0;
//                System.out.println("M0["+i+"]["+j+"] = "+M0[i]);
                a11[i][j] = h[i] * M1[i] / 2.0 + G1[i][j] / h[i] - M1[i] * h[i] / 6.0 - G1[i - 1][j] / h[i] + M1[i - 1] * h[i] / 6.0;
//                System.out.println("M1["+i+"]["+j+"] = "+M1[i]);
                a12[i][j] = h[i] * M2[i] / 2.0 + G2[i][j] / h[i] - M2[i] * h[i] / 6.0 - G2[i - 1][j] / h[i] + M2[i - 1] * h[i] / 6.0;
//                System.out.println("M2["+i+"]["+j+"] = "+M2[i]);
                a13[i][j] = h[i] * M3[i] / 2.0 + G3[i][j] / h[i] - M3[i] * h[i] / 6.0 - G3[i - 1][j] / h[i] + M3[i - 1] * h[i] / 6.0;
//                System.out.println("M3["+i+"]["+j+"] = "+M3[i]);

                a20[i][j] = M0[i];
                a21[i][j] = M1[i];
                a22[i][j] = M2[i];
                a23[i][j] = M3[i];

                a30[i][j] = (M0[i] - M0[i - 1]) / h[i];
                a31[i][j] = (M1[i] - M1[i - 1]) / h[i];
                a32[i][j] = (M2[i] - M2[i - 1]) / h[i];
                a33[i][j] = (M3[i] - M3[i - 1]) / h[i];
            }
        }
        for (int j=0;j<M+1;j++)
        {

            a00[0][j] = G0[0][j];//h[i] * M0[i] / 6.0 + G0[i][j] - M0[i] * h[i] * h[i] / 6.0;
            a01[0][j] = G1[0][j];//h[i] * M1[i] / 6.0 + G1[i][j] - M1[i] * h[i] * h[i] / 6.0;
            a02[0][j] = G2[0][j];//h[i] * M2[i] / 6.0 + G2[i][j] - M2[i] * h[i] * h[i] / 6.0;
            a03[0][j] = G3[0][j];//h[i] * M3[i] / 6.0 + G3[i][j] - M3[i] * h[i] * h[i] / 6.0;


            a20[0][j] = M0[0];
            a21[0][j] = M1[0];
            a22[0][j] = M2[0];
            a23[0][j] = M3[0];

//            a10[0][j] = 0.0;//Math.PI/2.0;
//            a11[0][j] = 0.0;
//            a12[0][j] = 0.0;
//            a13[0][j] = 0.0;
//
//            a30[0][j] = 0.0;//-Math.PI*Math.PI*Math.PI/8.0;
//            a31[0][j] = 0.0;
//            a32[0][j] = 0.0;
//            a33[0][j] = 0.0;




            a10[0][j] = a0y10.apply(uzelY[j]);
//            a10[0][j] = -1.0/epsilon;
            a11[0][j] = 0.0;
            a12[0][j] = 0.0;
            a13[0][j] = 0.0;

            a30[0][j] = a0y30.apply(uzelY[j]);
//            a30[0][j] = -1.0/(epsilon*epsilon*epsilon);
            a31[0][j] = 0.0;
            a32[0][j] = 0.0;
            a33[0][j] = 0.0;

        }


//        List<Double> norm = new ArrayList<>();
        findU(L,K,N,M,uzelX,uzelY,x,y,a00,a01,a02,a03,a10,a11,a12,a13,a20,a21,a22,a23,a30,a31,a32,a33,u,hGust,tGust);
//        x[0] = 0.;
//        for (int k = 1; k < L + 1; k++) {
//            x[k] = x[k - 1] + hSetkaShishkinaGust[k];
//        }
//        y[0] = 0.;
//        for (int m = 1; m < K + 1; m++) {
//            y[m] = y[m - 1] + tSetkaShishkinaGust[m];
//        }
        pogreshnost2(res2, L, K, x, y, u, biFunction, norm);


//
//        for(int i=0;i<N+1;i++){
//            double finalI = Double.valueOf(i);
//            Function<Double, Double> function = x -> biFunction.apply(x, finalI);
////            cubicalSpline(N, function, epsilon, norm);
//        }
//        System.out.println("Максимальная норма погрешности " + Collections.max(norm));
//        System.out.println("Минимальная норма погрешности " + Collections.min(norm));

        return u;
    }

    public static void findU (int L, int K, int N, int M, double[] uzelX, double[] uzelY, double[] x, double[] y,
                              double[][] a00, double[][] a01, double[][] a02, double[][] a03, double[][] a10, double[][] a11, double[][] a12,
                              double[][] a13, double[][] a20, double[][] a21, double[][] a22, double[][] a23, double[][] a30, double[][] a31,
                              double[][] a32, double[][] a33, double[][] U, double[] hGust, double[] tGust)
    {

        x[0] = 0.;
        for (int k = 1; k < L + 1; k++) {
            x[k] = x[k - 1] + hGust[k];
        }
        y[0] = 0.;
        for (int m = 1; m < K + 1; m++) {
            y[m] = y[m - 1] + tGust[m];
        }

        for (int i = 0; i < N ; i++) {
            for (int j = 0; j < M ; j++) {
                for (int k = 0; k < L +1; k++) {
                    for (int m = 0; m < K +1; m++) {
                        if ((x[k] >= uzelX[i]) && (x[k] < uzelX[i + 1]) && (y[m] >= uzelY[j]) && (y[m] < uzelY[j + 1])) {
//                                S[k] = (uzel[i + 1] - x[k]) * (uzel[i + 1] - x[k]) * (uzel[i + 1] - x[k]) * c[i] / (6. * h[i + 1]) + (x[k] - uzel[i]) * (x[k] - uzel[i]) * (x[k] - uzel[i]) * c[i + 1] / (6. * h[i + 1]) + (f[i + 1] / h[i + 1] - c[i + 1] * h[i + 1] / 6.) * (x[k] - uzel[i]) + (f[i] / h[i + 1] - c[i] * h[i + 1] / 6.) * (uzel[i + 1] - x[k]);
                            U[k][m] = a00[i][j] + a01[i][j] * (y[m] - uzelY[j]) + a02[i][j] * (y[m] - uzelY[j]) * (y[m] - uzelY[j]) / 2.0 + a03[i][j] * (y[m] - uzelY[j]) * (y[m] - uzelY[j]) * (y[m] - uzelY[j]) / 6.0 +
                                    +a10[i][j] * (x[k] - uzelX[i]) + a11[i][j] * (x[k] - uzelX[i]) * (y[m] - uzelY[j]) + a12[i][j] * (x[k] - uzelX[i]) * (y[m] - uzelY[j]) * (y[m] - uzelY[j]) / 2.0 +
                                    +a13[i][j] * (x[k] - uzelX[i]) * (y[m] - uzelY[j]) * (y[m] - uzelY[j]) * (y[m] - uzelY[j]) / 6.0 + a20[i][j] * (x[k] - uzelX[i]) * (x[k] - uzelX[i]) / 2.0 +
                                    +a21[i][j] * (x[k] - uzelX[i]) * (x[k] - uzelX[i]) * (y[m] - uzelY[j]) / 2.0 + a22[i][j] * (x[k] - uzelX[i]) * (x[k] - uzelX[i]) * (y[m] - uzelY[j]) * (y[m] - uzelY[j]) / 4.0 +
                                    +a23[i][j] * (x[k] - uzelX[i]) * (x[k] - uzelX[i]) * (y[m] - uzelY[j]) * (y[m] - uzelY[j]) * (y[m] - uzelY[j]) / 12.0 + a30[i][j] * (x[k] - uzelX[i]) * (x[k] - uzelX[i]) * (x[k] - uzelX[i]) / 6.0 +
                                    +a31[i][j] * (x[k] - uzelX[i]) * (x[k] - uzelX[i]) * (x[k] - uzelX[i]) * (y[m] - uzelY[j]) / 6.0 + a32[i][j] * (x[k] - uzelX[i]) * (x[k] - uzelX[i]) * (x[k] - uzelX[i]) * (y[m] - uzelY[j]) * (y[m] - uzelY[j]) / 12.0 +
                                    +a33[i][j] * (y[m] - (x[k] - uzelX[i]) * (x[k] - uzelX[i]) * (x[k] - uzelX[i]) * uzelY[j]) * (y[m] - uzelY[j]) * (y[m] - uzelY[j]) / 36.0;
//                            System.out.println("u["+k+","+m+"] = "+U[k][m]);
                        }
                    }
                }
            }
//            System.out.println();

//            S[L] = (uzel[N] - x[L])*(uzel[N] - x[L])*(uzel[N] - x[L])*c[N]/(6.*h[N]) + (x[L] - uzel[N-1])*(x[L] - uzel[N-1])*(x[L] - uzel[N-1])*c[N]/(6.*h[N]) + (f[N]/h[N] - c[N]*h[N]/6.)*(x[L] - uzel[N-1]) + (f[N-1]/h[N] - c[N-1]*h[N]/6.)*(uzel[N] - x[L]);

//            U[L][j] =
        }
        for (int j=0;j<M;j++)
        {
            for (int m=0;m<K+1;m++)
            {
                if ( (y[m] >= uzelY[j]) && (y[m] < uzelY[j + 1]) )
                {
                    U[L][m] = a00[N][j] + a01[N][j] * (y[m] - uzelY[j]) + a02[N][j] * (y[m] - uzelY[j]) * (y[m] - uzelY[j]) / 2.0 + a03[N][j] * (y[m] - uzelY[j]) * (y[m] - uzelY[j]) * (y[m] - uzelY[j]) / 6.0 +
                            +a10[N][j] * (x[L] - uzelX[N]) + a11[N][j] * (x[L] - uzelX[N]) * (y[m] - uzelY[j]) + a12[N][j] * (x[L] - uzelX[N]) * (y[m] - uzelY[j]) * (y[m] - uzelY[j]) / 2.0 +
                            +a13[N][j] * (x[L] - uzelX[N]) * (y[m] - uzelY[j]) * (y[m] - uzelY[j]) * (y[m] - uzelY[j]) / 6.0 + a20[N][j] * (x[L] - uzelX[N]) * (x[L] - uzelX[N]) / 2.0 +
                            +a21[N][j] * (x[L] - uzelX[N]) * (x[L] - uzelX[N]) * (y[m] - uzelY[j]) / 2.0 + a22[N][j] * (x[L] - uzelX[N]) * (x[L] - uzelX[N]) * (y[m] - uzelY[j]) * (y[m] - uzelY[j]) / 4.0 +
                            +a23[N][j] * (x[L] - uzelX[N]) * (x[L] - uzelX[N]) * (y[m] - uzelY[j]) * (y[m] - uzelY[j]) * (y[m] - uzelY[j]) / 12.0 + a30[N][j] * (x[L] - uzelX[N]) * (x[L] - uzelX[N]) * (x[L] - uzelX[N]) / 6.0 +
                            +a31[N][j] * (x[L] - uzelX[N]) * (x[L] - uzelX[N]) * (x[L] - uzelX[N]) * (y[m] - uzelY[j]) / 6.0 + a32[N][j] * (x[L] - uzelX[N]) * (x[L] - uzelX[N]) * (x[L] - uzelX[N]) * (y[m] - uzelY[j]) * (y[m] - uzelY[j]) / 12.0 +
                            +a33[N][j] * (y[m] - (x[L] - uzelX[N]) * (x[L] - uzelX[N]) * (x[L] - uzelX[N]) * uzelY[j]) * (y[m] - uzelY[j]) * (y[m] - uzelY[j]) / 36.0;
                }
            }
        }
        for (int i=0;i<N;i++)
        {
            for (int k=0;k<L+1;k++)
            {
                if ( (x[k] >= uzelX[i]) && (x[k] < uzelX[i + 1]) )
                {
                    U[k][K] = a00[i][M] + a01[i][M] * (y[K] - uzelY[M]) + a02[i][M] * (y[K] - uzelY[M]) * (y[K] - uzelY[M]) / 2.0 + a03[i][M] * (y[K] - uzelY[M]) * (y[K] - uzelY[M]) * (y[K] - uzelY[M]) / 6.0 +
                            +a10[i][M] * (x[k] - uzelX[i]) + a11[i][M] * (x[k] - uzelX[i]) * (y[K] - uzelY[M]) + a12[i][M] * (x[k] - uzelX[i]) * (y[K] - uzelY[M]) * (y[K] - uzelY[M]) / 2.0 +
                            +a13[i][M] * (x[k] - uzelX[i]) * (y[K] - uzelY[M]) * (y[K] - uzelY[M]) * (y[K] - uzelY[M]) / 6.0 + a20[i][M] * (x[k] - uzelX[i]) * (x[k] - uzelX[i]) / 2.0 +
                            +a21[i][M] * (x[k] - uzelX[i]) * (x[k] - uzelX[i]) * (y[K] - uzelY[M]) / 2.0 + a22[i][M] * (x[k] - uzelX[i]) * (x[k] - uzelX[i]) * (y[K] - uzelY[M]) * (y[K] - uzelY[M]) / 4.0 +
                            +a23[i][M] * (x[k] - uzelX[i]) * (x[k] - uzelX[i]) * (y[K] - uzelY[M]) * (y[K] - uzelY[M]) * (y[K] - uzelY[M]) / 12.0 + a30[i][M] * (x[k] - uzelX[i]) * (x[k] - uzelX[i]) * (x[k] - uzelX[i]) / 6.0 +
                            +a31[i][M] * (x[k] - uzelX[i]) * (x[k] - uzelX[i]) * (x[k] - uzelX[i]) * (y[K] - uzelY[M]) / 6.0 + a32[i][M] * (x[k] - uzelX[i]) * (x[k] - uzelX[i]) * (x[k] - uzelX[i]) * (y[K] - uzelY[M]) * (y[K] - uzelY[M]) / 12.0 +
                            +a33[i][M] * (y[K] - (x[k] - uzelX[i]) * (x[k] - uzelX[i]) * (x[k] - uzelX[i]) * uzelY[M]) * (y[K] - uzelY[M]) * (y[K] - uzelY[M]) / 36.0;
                }
            }
        }
        U[L][K] = a00[N][M] + a01[N][M] * (y[K] - uzelY[M]) + a02[N][M] * (y[K] - uzelY[M]) * (y[K] - uzelY[M]) / 2.0 + a03[N][M] * (y[K] - uzelY[M]) * (y[K] - uzelY[M]) * (y[K] - uzelY[M]) / 6.0 +
                +a10[N][M] * (x[L] - uzelX[N]) + a11[N][M] * (x[L] - uzelX[N]) * (y[K] - uzelY[M]) + a12[N][M] * (x[L] - uzelX[N]) * (y[K] - uzelY[M]) * (y[K] - uzelY[M]) / 2.0 +
                +a13[N][M] * (x[L] - uzelX[N]) * (y[K] - uzelY[M]) * (y[K] - uzelY[M]) * (y[K] - uzelY[M]) / 6.0 + a20[N][M] * (x[L] - uzelX[N]) * (x[L] - uzelX[N]) / 2.0 +
                +a21[N][M] * (x[L] - uzelX[N]) * (x[L] - uzelX[N]) * (y[K] - uzelY[M]) / 2.0 + a22[N][M] * (x[L] - uzelX[N]) * (x[L] - uzelX[N]) * (y[K] - uzelY[M]) * (y[K] - uzelY[M]) / 4.0 +
                +a23[N][M] * (x[L] - uzelX[N]) * (x[L] - uzelX[N]) * (y[K] - uzelY[M]) * (y[K] - uzelY[M]) * (y[K] - uzelY[M]) / 12.0 + a30[N][M] * (x[L] - uzelX[N]) * (x[L] - uzelX[N]) * (x[L] - uzelX[N]) / 6.0 +
                +a31[N][M] * (x[L] - uzelX[N]) * (x[L] - uzelX[N]) * (x[L] - uzelX[N]) * (y[K] - uzelY[M]) / 6.0 + a32[N][M] * (x[L] - uzelX[N]) * (x[L] - uzelX[N]) * (x[L] - uzelX[N]) * (y[K] - uzelY[M]) * (y[K] - uzelY[M]) / 12.0 +
                +a33[N][M] * (y[K] - (x[L] - uzelX[N]) * (x[L] - uzelX[N]) * (x[L] - uzelX[N]) * uzelY[M]) * (y[K] - uzelY[M]) * (y[K] - uzelY[M]) / 36.0;

    }

//        public static void main(String[] args){
//            double epsilon = 1.e-4;
//            System.out.println("Epsilon = "+epsilon);
//            double a = 0.;
//            for (int i=16; i<=4096; i=i*2){
//                List<Double> norm = new ArrayList<>();
//                cubicalSpline(i, x->Math.cos(Math.PI*x/2.) + Math.exp(-x/epsilon), epsilon, norm, eps->-Math.PI*Math.PI/4.+1.0/(eps*eps), eps->Math.exp(-1./eps)/(eps*eps));
//                double four = a/norm.get(0);
//                a = norm.get(0);
//                System.out.println("Порядок точности log2 (||"+i/2.+"||/||"+i+"||) = "+Math.log10(four)/Math.log10(2.));
//                System.out.println("_______________________________________________________________");
//                System.out.println("||"+i+"|| = "+ norm);
//
//            }
//        }

    public static void main(String[] args){
        double epsilon = 1.e-03; //22:49;
//        double epsilonScannerDouble;
        String equalation = "d_first_derivation_plus_e-y";
        String setkaX = "r";
        String setkaY = "r";


//        System.out.print("Enter the epsilon = ");
//        Scanner epsilonScanner = new Scanner(System.in);
//        epsilonScannerDouble = epsilonScanner.nextDouble();
//        System.out.println();
//        System.out.print("Mesh for X = ");
//        Scanner mesh1Scanner = new Scanner(System.in);
//        String setkaX = mesh1Scanner.next();
//        System.out.println("Mesh for Y = ");
//        Scanner mesh2Scanner = new Scanner(System.in);
//        double epsilon = epsilonScannerDouble;
//        String setkaY = mesh2Scanner.next();


        System.out.println(setkaX);
        System.out.println(setkaY);
        System.out.println("Epsilon = "+epsilon);
        double a = 0.;
        for (int i=64; i<=512; i=i*2){
            List<Double> norm = new ArrayList<>();
            BiFunction<Double, Double, Double> equal = null;
            Function<Double, Double> secondDerivationY0 = null;
            Function<Double, Double> secondDerivationYN= null;
            Function<Double, Double> betta1Y= null;
            Function<Double, Double> nullSecondDerivationX0= null;
            Function<Double, Double> nullSecondDerivationXN= null;
            Function<Double, Double> nullBetta1X= null;
            Function<Double, Double> firstSecondDerivationX0= null;
            Function<Double, Double> firstSecondDerivationXN= null;
            Function<Double, Double> firstBetta1X= null;
            Function<Double, Double> secondSecondDerivationX0= null;
            Function<Double, Double> secondSecondDerivationXN= null;
            Function<Double, Double> secondBetta1X= null;
            Function<Double, Double> thirdSecondDerivationX0= null;
            Function<Double, Double> thirdSecondDerivationXN= null;
            Function<Double, Double> thirdBetta1X= null;
            Function<Double, Double> G1i0= null;
            Function<Double, Double> G3i0= null;
            Function<Double, Double> a0y10= null;
            Function<Double, Double> a0y30= null;







            if (equalation == "a") {
                equal = (x, y) -> Math.cos(Math.PI * x / 2.) + Math.exp(-y / epsilon);

                secondDerivationY0 = x -> 1.0 / (epsilon * epsilon);
                secondDerivationYN = x -> Math.exp(-1.0 / epsilon) / (epsilon * epsilon);
                betta1Y = x -> 1.0 / (epsilon * epsilon);
                nullSecondDerivationX0 = y -> -Math.PI * Math.PI / (4.0);
                nullSecondDerivationXN = y -> 0.0;
                nullBetta1X = y -> -Math.PI * Math.PI / (4.0);
                firstSecondDerivationX0 = y -> 0.0;
                firstSecondDerivationXN = y -> 0.0;
                firstBetta1X = y -> 0.0;
                secondSecondDerivationX0 = y -> 0.0;
                secondSecondDerivationXN = y -> 0.0;
                secondBetta1X = y -> 0.0;
                thirdSecondDerivationX0 = y -> 0.0;
                thirdSecondDerivationXN = y -> 0.0;
                thirdBetta1X = y -> 0.0;
                G1i0 = x -> -1.0 / epsilon;
                G3i0 = x -> -1.0 / (epsilon * epsilon * epsilon);
                a0y10 = y -> 0.0;
                a0y30 = y -> 0.0;


//                equal = (x, y) -> Math.cos(Math.PI * x / 2.) + Math.exp(-y / epsilon);
//
//                secondDerivationY0 = x -> 1.0 / (epsilon);
//                secondDerivationYN = x -> -Math.exp(-1.0 / epsilon) / (epsilon);
//                betta1Y = x -> 1.0 / (epsilon * epsilon);
//                nullSecondDerivationX0 = y -> 0.0;
//                nullSecondDerivationXN = y -> -Math.PI/2.0;
//                nullBetta1X = y -> -Math.PI * Math.PI / (4.0);
//                firstSecondDerivationX0 = y -> 0.0;
//                firstSecondDerivationXN = y -> 0.0;
//                firstBetta1X = y -> 0.0;
//                secondSecondDerivationX0 = y -> 0.0;
//                secondSecondDerivationXN = y -> 0.0;
//                secondBetta1X = y -> 0.0;
//                thirdSecondDerivationX0 = y -> 0.0;
//                thirdSecondDerivationXN = y -> 0.0;
//                thirdBetta1X = y -> 0.0;
//                G1i0 = x -> -1.0 / epsilon;
//                G3i0 = x -> -1.0 / (epsilon * epsilon * epsilon);
//                a0y10 = y -> 0.0;
//                a0y30 = y -> 0.0;
            }

            if (equalation == "b") {
                equal = (x, y) -> Math.exp(-x / epsilon) + Math.exp(-2.0 * y / epsilon) + Math.cos(Math.PI * x / 2.0) * Math.exp(-2.0 * y);

                secondDerivationY0 = x -> 4.0 / (epsilon * epsilon* epsilon) + 4.0 * Math.cos(Math.PI * x / 2.0);
                secondDerivationYN = x -> 4.0 / (epsilon * epsilon* epsilon) * Math.exp(-2.0 / epsilon) + 4.0 * Math.cos(Math.PI * x / 2.0) * Math.exp(-2.0);
                betta1Y = x -> 4.0 / (epsilon * epsilon* epsilon) + 4.0 * Math.cos(Math.PI * x / 2.0);

                nullSecondDerivationX0 = y -> 1.0 / (epsilon * epsilon) - Math.PI * Math.PI / 4.0 * Math.exp(-2.0 * y);
                nullSecondDerivationXN = y -> 1.0 / (epsilon * epsilon) * Math.exp(-1.0 / epsilon);
                nullBetta1X = y -> 1.0 / (epsilon * epsilon) - Math.PI * Math.PI / 4.0 * Math.exp(-2.0 * y);

                firstSecondDerivationX0 = y -> Math.PI * Math.PI / 2.0 * Math.exp(-2.0 * y);
                firstSecondDerivationXN = y -> 0.0;
                firstBetta1X = y -> Math.PI * Math.PI / 2.0 * Math.exp(-2.0 * y);

                secondSecondDerivationX0 = y -> -Math.PI * Math.PI * Math.exp(-2.0 * y);
                secondSecondDerivationXN = y -> 0.0;
                secondBetta1X = y -> -Math.PI * Math.PI * Math.exp(-2.0 * y);

                thirdSecondDerivationX0 = y -> 2.0 * Math.PI * Math.PI * Math.exp(-2.0 * y);
                thirdSecondDerivationXN = y -> 0.0;
                thirdBetta1X = y -> 2.0 * Math.PI * Math.PI * Math.exp(-2.0 * y);

                G1i0 = x -> -2.0 / epsilon - 2.0 * Math.cos(Math.PI * x / 2.0);
                G3i0 = x -> -8.0 / (epsilon * epsilon * epsilon) - 8.0 * Math.cos(Math.PI * x / 2.0);
                a0y10 = y -> -1.0 / epsilon;
                a0y30 = y -> -1.0 / (epsilon * epsilon * epsilon);
            }

            if (equalation == "c") {
                equal = (x, y) -> Math.exp(-x/epsilon) + Math.exp(-y/epsilon) + Math.cos(Math.PI * x / 2.0) ;

                secondDerivationY0 = x -> 1.0 / (epsilon * epsilon);
                secondDerivationYN = x -> 1.0 / (epsilon * epsilon) * Math.exp(-1.0 / epsilon);
                betta1Y = x -> 1.0 / (epsilon * epsilon);

                nullSecondDerivationX0 = y -> 1.0 / (epsilon * epsilon) - Math.PI*Math.PI/4.0;
                nullSecondDerivationXN = y -> 1.0 / (epsilon * epsilon) * Math.exp(-1.0 / epsilon);
                nullBetta1X = y -> 1.0 / (epsilon * epsilon) - Math.PI*Math.PI/4.0;

                firstSecondDerivationX0 = y -> 0.0;
                firstSecondDerivationXN = y -> 0.0;
                firstBetta1X = y -> 0.0;

                secondSecondDerivationX0 = y -> 0.0;
                secondSecondDerivationXN = y -> 0.0;
                secondBetta1X = y -> 0.0;

                thirdSecondDerivationX0 = y -> 0.0;
                thirdSecondDerivationXN = y -> 0.0;
                thirdBetta1X = y -> 0.0;

                G1i0 = x -> -1.0 / epsilon;
                G3i0 = x -> -1.0 / (epsilon * epsilon * epsilon) ;
                a0y10 = y -> -1.0 / epsilon;
                a0y30 = y -> -1.0 / (epsilon * epsilon * epsilon);
            }

            if (equalation == "c_exp") {
                equal = (x, y) -> Math.exp(-x/epsilon) + Math.exp(-y/epsilon) + 5.0*Math.cos(Math.PI * x / 2.0) ;

                secondDerivationY0 = x -> 1.0 / (epsilon * epsilon);
                secondDerivationYN = x -> 1.0 / (epsilon * epsilon) * Math.exp(-1.0 / epsilon);
                betta1Y = x -> 1.0 / (epsilon * epsilon);

                nullSecondDerivationX0 = y -> 1.0 / (epsilon * epsilon) - 5*Math.PI*Math.PI/4.0;
                nullSecondDerivationXN = y -> 1.0 / (epsilon * epsilon) * Math.exp(-1.0 / epsilon);
                nullBetta1X = y -> 1.0 / (epsilon * epsilon) - 5*Math.PI*Math.PI/4.0;

                firstSecondDerivationX0 = y -> 0.0;
                firstSecondDerivationXN = y -> 0.0;
                firstBetta1X = y -> 0.0;

                secondSecondDerivationX0 = y -> 0.0;
                secondSecondDerivationXN = y -> 0.0;
                secondBetta1X = y -> 0.0;

                thirdSecondDerivationX0 = y -> 0.0;
                thirdSecondDerivationXN = y -> 0.0;
                thirdBetta1X = y -> 0.0;

                G1i0 = x -> -1.0 / epsilon;
                G3i0 = x -> -1.0 / (epsilon * epsilon * epsilon) ;
                a0y10 = y -> -1.0 / epsilon;
                a0y30 = y -> -1.0 / (epsilon * epsilon * epsilon);
            }



            if (equalation == "d") {
                equal = (x, y) -> Math.exp(-x/epsilon) + Math.exp(-2.0*y/epsilon) + Math.cos(Math.PI * x / 2.0) ;

                secondDerivationY0 = x -> 4.0 / (epsilon * epsilon);
                secondDerivationYN = x -> 4.0 / (epsilon * epsilon) * Math.exp(-2.0 / epsilon);
                betta1Y = x -> 4.0 / (epsilon * epsilon);

                nullSecondDerivationX0 = y -> 1.0 / (epsilon * epsilon) - Math.PI*Math.PI/4.0;
                nullSecondDerivationXN = y -> 1.0 / (epsilon * epsilon) * Math.exp(-1.0 / epsilon);
                nullBetta1X = y -> 1.0 / (epsilon * epsilon) - Math.PI*Math.PI/4.0;

                firstSecondDerivationX0 = y -> 0.0;
                firstSecondDerivationXN = y -> 0.0;
                firstBetta1X = y -> 0.0;

                secondSecondDerivationX0 = y -> 0.0;
                secondSecondDerivationXN = y -> 0.0;
                secondBetta1X = y -> 0.0;

                thirdSecondDerivationX0 = y -> 0.0;
                thirdSecondDerivationXN = y -> 0.0;
                thirdBetta1X = y -> 0.0;

                G1i0 = x -> -2.0 / epsilon;
                G3i0 = x -> -8.0 / (epsilon * epsilon * epsilon) ;
                a0y10 = y -> -1.0 / epsilon;
                a0y30 = y -> -1.0 / (epsilon * epsilon * epsilon);
            }


            if (equalation == "d_first_derivation") {
                equal = (x, y) -> Math.exp(-x/epsilon) + Math.exp(-2.0*y/epsilon) + Math.cos(Math.PI * x / 2.0) ;

                secondDerivationY0 = x -> -2.0 / (epsilon);
                secondDerivationYN = x -> -2.0 / (epsilon) * Math.exp(-2.0 / epsilon);
                betta1Y = x -> 4.0 / (epsilon * epsilon);

                nullSecondDerivationX0 = y -> -1.0 / (epsilon);
                nullSecondDerivationXN = y -> -Math.PI/2.0 - 1.0 / (epsilon) * Math.exp(-1.0 / epsilon);
                nullBetta1X = y -> 1.0 / (epsilon * epsilon) - Math.PI*Math.PI/4.0;

                firstSecondDerivationX0 = y -> 0.0;
                firstSecondDerivationXN = y -> 0.0;
                firstBetta1X = y -> 0.0;

                secondSecondDerivationX0 = y -> 0.0;
                secondSecondDerivationXN = y -> 0.0;
                secondBetta1X = y -> 0.0;

                thirdSecondDerivationX0 = y -> 0.0;
                thirdSecondDerivationXN = y -> 0.0;
                thirdBetta1X = y -> 0.0;

                G1i0 = x -> -2.0 / epsilon;
                G3i0 = x -> -8.0 / (epsilon * epsilon * epsilon) ;
                a0y10 = y -> -1.0 / epsilon;
                a0y30 = y -> -1.0 / (epsilon * epsilon * epsilon);
            }

            if (equalation == "bs") {
                equal = (x, y) -> Math.exp(-x / epsilon) + Math.exp(-2.0 * y / epsilon) + Math.cos(Math.PI * x / 2.0) * Math.exp(-y);

                secondDerivationY0 = x -> 4.0 / (epsilon * epsilon* epsilon) + Math.cos(Math.PI * x / 2.0);
                secondDerivationYN = x -> 4.0 / (epsilon * epsilon* epsilon) * Math.exp(-2.0 / epsilon) + Math.cos(Math.PI * x / 2.0) * Math.exp(-1.0);
                betta1Y = x -> 4.0 / (epsilon * epsilon* epsilon) + Math.cos(Math.PI * x / 2.0);

                nullSecondDerivationX0 = y -> 1.0 / (epsilon * epsilon) - Math.PI * Math.PI / 4.0 * Math.exp(-y);
                nullSecondDerivationXN = y -> 1.0 / (epsilon * epsilon) * Math.exp(-1.0 / epsilon);
                nullBetta1X = y -> 1.0 / (epsilon * epsilon) - Math.PI * Math.PI / 4.0 * Math.exp(-y);

                firstSecondDerivationX0 = y -> Math.PI * Math.PI / 4.0 * Math.exp(-y);
                firstSecondDerivationXN = y -> 0.0;
                firstBetta1X = y -> Math.PI * Math.PI / 4.0 * Math.exp(-y);

                secondSecondDerivationX0 = y -> -Math.PI * Math.PI/4.0 * Math.exp(-y);
                secondSecondDerivationXN = y -> 0.0;
                secondBetta1X = y -> -Math.PI * Math.PI/4.0 * Math.exp(-y);

                thirdSecondDerivationX0 = y -> Math.PI * Math.PI/4.0 * Math.exp(-y);
                thirdSecondDerivationXN = y -> 0.0;
                thirdBetta1X = y ->  Math.PI * Math.PI/4.0 * Math.exp(-y);

                G1i0 = x -> -2.0 / epsilon - Math.cos(Math.PI * x / 2.0);
                G3i0 = x -> -8.0 / (epsilon * epsilon * epsilon) - Math.cos(Math.PI * x / 2.0);
                a0y10 = y -> -1.0 / epsilon;
                a0y30 = y -> -1.0 / (epsilon * epsilon * epsilon);
            }

            if (equalation == "bsin") {
                equal = (x, y) -> Math.exp(-x / epsilon) + Math.exp(-2.0 * y / epsilon) + Math.sin(Math.PI * x / 2.0) * Math.exp(-y);

                secondDerivationY0 = x -> 4.0 / (epsilon * epsilon) + Math.sin(Math.PI * x / 2.0);
                secondDerivationYN = x -> 4.0 / (epsilon * epsilon) * Math.exp(-2.0 / epsilon) + Math.sin(Math.PI * x / 2.0) * Math.exp(-1.0);
                betta1Y = x -> 4.0 / (epsilon * epsilon) + Math.sin(Math.PI * x / 2.0);

                nullSecondDerivationX0 = y -> 1.0 / (epsilon * epsilon) ;
                nullSecondDerivationXN = y -> 1.0 / (epsilon * epsilon)*Math.exp(-1.0/epsilon) - Math.PI * Math.PI / 4.0 * Math.exp(-y);
                nullBetta1X = y -> 1.0 / (epsilon * epsilon) ;

                firstSecondDerivationX0 = y -> 0.0;
                firstSecondDerivationXN = y -> Math.PI * Math.PI / 4.0 * Math.exp(-y);
                firstBetta1X = y -> 0.0;

                secondSecondDerivationX0 = y -> 0.0;
                secondSecondDerivationXN = y -> -Math.PI * Math.PI/4.0 * Math.exp(-y);
                secondBetta1X = y -> 0.0;

                thirdSecondDerivationX0 = y -> 0.0;
                thirdSecondDerivationXN = y -> Math.PI * Math.PI/4.0 * Math.exp(-y);
                thirdBetta1X = y ->  0.0;

                G1i0 = x -> -2.0 / epsilon - Math.sin(Math.PI * x / 2.0);
                G3i0 = x -> -8.0 / (epsilon * epsilon * epsilon) - Math.sin(Math.PI * x / 2.0);
                a0y10 = y -> -1.0 / epsilon + Math.PI/2.0*Math.exp(-y);
                a0y30 = y -> -1.0 / (epsilon * epsilon * epsilon) - Math.PI*Math.PI*Math.PI/8.0*Math.exp(-y);
            }


            if (equalation == "bsinexperimental") {
                equal = (x, y) -> Math.exp(-x / epsilon) + Math.exp(-2.0 * y / epsilon) + Math.sin(Math.PI * x / 2.0) * Math.exp(-y);

                secondDerivationY0 = x -> 4.0 / (epsilon * epsilon) + Math.sin(Math.PI * x / 2.0);
                secondDerivationYN = x -> 4.0 / (epsilon * epsilon) * Math.exp(-2.0 / epsilon) + Math.sin(Math.PI * x / 2.0) * Math.exp(-1.0);
                betta1Y = x -> 4.0 / (epsilon * epsilon) + Math.sin(Math.PI * x / 2.0);

                nullSecondDerivationX0 = y -> 1.0 / (epsilon * epsilon) ;
                nullSecondDerivationXN = y -> 1.0 / (epsilon * epsilon)*Math.exp(-1.0/epsilon) - Math.PI * Math.PI / 4.0 * Math.exp(-y);
                nullBetta1X = y -> 1.0 / (epsilon * epsilon) ;

                firstSecondDerivationX0 = y -> 0.0;
                firstSecondDerivationXN = y -> Math.PI * Math.PI / 4.0 * Math.exp(-y);
                firstBetta1X = y -> 0.0;

                secondSecondDerivationX0 = y -> 0.0;
                secondSecondDerivationXN = y -> -Math.PI * Math.PI/4.0 * Math.exp(-y);
                secondBetta1X = y -> 0.0;

                thirdSecondDerivationX0 = y -> 0.0;
                thirdSecondDerivationXN = y -> Math.PI * Math.PI/4.0 * Math.exp(-y);
                thirdBetta1X = y ->  0.0;

                G1i0 = x -> -2.0 / epsilon - Math.sin(Math.PI * x / 2.0);
                G3i0 = x -> -8.0 / (epsilon * epsilon * epsilon) - Math.sin(Math.PI * x / 2.0);
                a0y10 = y -> -1.0 / epsilon + Math.PI/2.0*Math.exp(-y);
                a0y30 = y -> -1.0 / (epsilon * epsilon * epsilon) - Math.PI*Math.PI*Math.PI/8.0*Math.exp(-y);
            }

            if (equalation == "d_first_derivation_plus_e-y") {
                equal = (x, y) -> Math.exp(-x/epsilon) + Math.exp(-2.0*y/epsilon) + Math.cos(Math.PI * x / 2.0) + Math.exp(y) ;

                secondDerivationY0 = x -> -2.0 / (epsilon) + 1.0;
                secondDerivationYN = x -> -2.0 / (epsilon) * Math.exp(-2.0 / epsilon) + Math.exp(1.0);
                betta1Y = x -> 4.0 / (epsilon * epsilon) + 1.0;

                nullSecondDerivationX0 = y -> -1.0 / (epsilon);
                nullSecondDerivationXN = y -> -Math.PI/2.0 - 1.0 / (epsilon) * Math.exp(-1.0 / epsilon);
                nullBetta1X = y -> 1.0 / (epsilon * epsilon) - Math.PI*Math.PI/4.0;

                firstSecondDerivationX0 = y -> 0.0;
                firstSecondDerivationXN = y -> 0.0;
                firstBetta1X = y -> 0.0;

                secondSecondDerivationX0 = y -> 0.0;
                secondSecondDerivationXN = y -> 0.0;
                secondBetta1X = y -> 0.0;

                thirdSecondDerivationX0 = y -> 0.0;
                thirdSecondDerivationXN = y -> 0.0;
                thirdBetta1X = y -> 0.0;

                G1i0 = x -> -2.0 / epsilon + 1.0;
                G3i0 = x -> -8.0 / (epsilon * epsilon * epsilon) + 1.0;
                a0y10 = y -> -1.0 / epsilon;
                a0y30 = y -> -1.0 / (epsilon * epsilon * epsilon);
            }



            if (equalation == "d_plus_e-y") {
                equal = (x, y) -> Math.exp(-x/epsilon) + Math.exp(-2.0*y/epsilon) + Math.cos(Math.PI * x / 2.0) + Math.exp(y);

                secondDerivationY0 = x -> 4.0 / (epsilon * epsilon) + 1.0;
                secondDerivationYN = x -> 4.0 / (epsilon * epsilon) * Math.exp(-2.0 / epsilon) + Math.exp(1.0);
                betta1Y = x -> 4.0 / (epsilon * epsilon) + 1.0;

                nullSecondDerivationX0 = y -> 1.0 / (epsilon * epsilon) - Math.PI*Math.PI/4.0;
                nullSecondDerivationXN = y -> 1.0 / (epsilon * epsilon) * Math.exp(-1.0 / epsilon);
                nullBetta1X = y -> 1.0 / (epsilon * epsilon) - Math.PI*Math.PI/4.0;

                firstSecondDerivationX0 = y -> 0.0;
                firstSecondDerivationXN = y -> 0.0;
                firstBetta1X = y -> 0.0;

                secondSecondDerivationX0 = y -> 0.0;
                secondSecondDerivationXN = y -> 0.0;
                secondBetta1X = y -> 0.0;

                thirdSecondDerivationX0 = y -> 0.0;
                thirdSecondDerivationXN = y -> 0.0;
                thirdBetta1X = y -> 0.0;

                G1i0 = x -> -2.0 / epsilon + 1.0;
                G3i0 = x -> -8.0 / (epsilon * epsilon * epsilon) + 1.0;
                a0y10 = y -> -1.0 / epsilon;
                a0y30 = y -> -1.0 / (epsilon * epsilon * epsilon);
            }






//            Function<Double, Double> secondDerivationY0 = eps -> -1.0 / (eps);
//            Function<Double, Double> secondDerivationYN = eps -> -Math.exp(-1.0 / eps) / (eps);
//            Function<Double, Double> betta1Y = eps -> 1.0 / (eps * eps);
//            Function<Double, Double> secondDerivationX0 = eps -> 0.0;
//            Function<Double, Double> secondDerivationXN = eps -> -Math.PI / (2.0);
//            Function<Double, Double> betta1X = eps -> -Math.PI * Math.PI / (4.0);

//            cubicalSpline(i, x->Math.cos(Math.PI*x/2.) + Math.exp(-x/epsilon), epsilon, norm, eps->-Math.PI*Math.PI/4.+1.0/(eps*eps), eps->Math.exp(-1./eps)/(eps*eps), setkaX);//вторая произв
//            cubicalSpline(i, x->Math.cos(Math.PI*x/2.) + Math.exp(-x/epsilon), epsilon, norm, eps->-1.0/(eps), eps->-Math.PI/2.0-Math.exp(-1./eps)/(eps), setkaX);//первая произв
            bicubicalSpline(i,i, 5*i, 5*i, epsilon, equal, secondDerivationY0, secondDerivationYN, betta1Y, nullSecondDerivationX0, nullSecondDerivationXN, nullBetta1X, firstSecondDerivationX0, firstSecondDerivationXN, firstBetta1X, secondSecondDerivationX0, secondSecondDerivationXN, secondBetta1X, thirdSecondDerivationX0, thirdSecondDerivationXN, thirdBetta1X, norm, setkaX, setkaY, G1i0, G3i0, a0y10, a0y30);
            double four = a/norm.get(0);
            a = norm.get(0);
            System.out.println("Порядок точности log2 (||"+i/2.+"||/||"+i+"||) = "+Math.log10(four)/Math.log10(2.));
            System.out.println("_______________________________________________________________");
//            System.out.println("||"+i+"|| = "+ norm.get(0));
            System.out.println("||"+i+"|| = "+ norm.get(0)+"  ["+norm.get(1)+"]"+"["+norm.get(2)+"]");

        }
    }



    }
