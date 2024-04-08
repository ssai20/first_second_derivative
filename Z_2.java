public class Z_2 {

    public static double findUexp (int N, double epsilon) {
        int L = 4*(N-1);
        double[] U = new double[L+1];
        double[] proizvU = new double[L+1];
        double[] proizvPogranSloiU = new double[L+1];
        double hh = 1./L;
        double[] x = new double[L+1];

        x[0] = 0;
        for (int j=1;j<=L;j++){
            x[j] = x[j-1]+hh;
        }

        for (int i=0;i<=L;i++){
            U[i] = Math.cos(Math.PI*x[i]) + Math.exp(-x[i]/epsilon);
        }
        for (int i=1;i<=L;i++){
            proizvU[i] = (U[i] - U[i-1])/hh;
            proizvPogranSloiU[i] = (-1./epsilon)*Math.exp(-x[i]/epsilon)*(U[i] - U[i-1])/
                    (Math.exp(-x[i]/epsilon) - Math.exp(-x[i-1]/epsilon));
        }
        //вычисление первой производной на 5 точках
        for (int i=2;i<=L-2;i=i+4){
            proizvU[i-2] = (U[i] - U[i-1])/hh;
            proizvU[i-1] = (U[i] - U[i-1])/hh;
            proizvU[i] = (U[i] - U[i-1])/hh;
            proizvU[i+1] = (U[i] - U[i-1])/hh;
            proizvU[i+2] = (U[i] - U[i-1])/hh;
        }

        double[] norm = new double[L+1];

        for (int i=2;i<=L-2;i=i+4){
            norm[i-2] = Math.exp((x[i-1]-x[i-2])/epsilon) * (Math.cos(Math.PI*x[i]) - Math.cos(Math.PI*x[i-1]))/(1.-Math.exp(-hh/epsilon)) + epsilon*Math.PI*Math.sin(Math.PI*x[i-2]);
            norm[i-1] = Math.exp((x[i-1]-x[i-1])/epsilon) * (Math.cos(Math.PI*x[i]) - Math.cos(Math.PI*x[i-1]))/(1.-Math.exp(-hh/epsilon)) + epsilon*Math.PI*Math.sin(Math.PI*x[i-1]);
            norm[i] = Math.exp((x[i-1]-x[i])/epsilon) * (Math.cos(Math.PI*x[i]) - Math.cos(Math.PI*x[i-1]))/(1.-Math.exp(-hh/epsilon)) + epsilon*Math.PI*Math.sin(Math.PI*x[i]);
            norm[i+1] = Math.exp((x[i-1]-x[i+1])/epsilon) * (Math.cos(Math.PI*x[i]) - Math.cos(Math.PI*x[i-1]))/(1.-Math.exp(-hh/epsilon)) + epsilon*Math.PI*Math.sin(Math.PI*x[i+1]);
            norm[i+2] = Math.exp((x[i-1]-x[i+2])/epsilon) * (Math.cos(Math.PI*x[i]) - Math.cos(Math.PI*x[i-1]))/(1.-Math.exp(-hh/epsilon)) + epsilon*Math.PI*Math.sin(Math.PI*x[i+2]);
        }


        double maxNorm = 0;

        for (int i=0;i<=L;i++){
//            norm[i] = epsilon*Math.abs(proizvU[i]-(-Math.PI*Math.sin(Math.PI*x[i]) - Math.exp(-x[i]/epsilon)/epsilon));

            /*norm[i] = epsilon*Math.abs(proizvPogranSloiU[i] -
                    (-Math.PI*Math.sin(Math.PI*x[i]) - Math.exp(-x[i]/epsilon)/epsilon));*/
        }

        for(int i=0;i<=L;i++){
            if (maxNorm<norm[i]) maxNorm = norm[i];
        }
        return maxNorm;
    }

    public static double findUsqrt (int N, double epsilon) {
        int L = 5*N;
        double[] U = new double[L+1];
        double[] proizvU = new double[L+1];
        double[] proizvPogranSloiU = new double[L+1];
        double hh = 1./L;
        double[] x = new double[L+1];

        System.out.println("first differentional x+epsilon:");
        x[0] = 0;
        for (int j=1;j<=L;j++){
            x[j] = x[j-1]+hh;
        }

        for (int i=0;i<=L;i++){
            U[i] = Math.cos(Math.PI*x[i]) + Math.exp(-x[i]/epsilon);
        }
        for (int i=1;i<=L;i++){
            proizvU[i] = (U[i] - U[i-1])/hh;
            proizvPogranSloiU[i] = (U[i] - U[i-1])/
                    ((Math.sqrt(x[i] + epsilon) - Math.sqrt(x[i-1] + epsilon))*(2.*Math.sqrt(x[i]+epsilon)));
        }
        double maxNorm = 0;
        double[] norm = new double[L+1];
        for (int i=1;i<=L;i++){
            norm[i] = Math.sqrt(epsilon)*Math.abs(proizvU[i]-(-Math.PI*Math.sin(Math.PI*x[i]) - Math.exp(-x[i]/epsilon)/epsilon));
//            norm[i] = epsilon*Math.abs(proizvPogranSloiU[i] -
//                    (-Math.PI*Math.sin(Math.PI*x[i]) - Math.exp(-x[i]/epsilon)/epsilon));
        }

        for(int i=1;i<=L;i++){
            if (maxNorm<norm[i]) maxNorm = norm[i];
        }
        return maxNorm;
    }




    public static double findSecondUexp (int N, double epsilon) {
        int L = 4*(N-1);
//        System.out.println(L);
        double[] U = new double[L+1];
        double[] Fi = new double[L+1];
        double[] secondProizvU = new double[L+1];
        double[] secondProizvPogranSloiU = new double[L+1];
//        double h = 1/N;
        double hh = 1./L;
        double[] uzelX = new double[N+1];
        double[] x = new double[L+1];

//        uzelX[0] = 0;
//        for (int i=1;i<=N;i++){
//            uzelX[i] = uzelX[i-1]+h;
//        }
        x[0] = 0;
        for (int j=1;j<=L;j++){
            x[j] = x[j-1]+hh;
//            System.out.println(x[j]);
        }

        for (int i=0;i<=L;i++){
            U[i] = Math.cos(Math.PI*x[i]) + Math.exp(-x[i]/epsilon);
            Fi[i] = Math.exp(-x[i]/epsilon);
//            System.out.println(U[i]);
        }
//        for (int i=1;i<=L-1;i=i+4){
//            secondProizvU[i-1] = (U[i+1] - U[i-1])/(hh*hh);
//            secondProizvU[i] = (U[i+1] -2*U[i] + U[i-1])/(hh*hh);
//            secondProizvU[i+1] = (U[i+1] - U[i-1]) /(hh*hh);
//
////            System.out.println(proizvU[i]);
////            secondProizvPogranSloiU[i-1] = (1./(epsilon*epsilon))*Math.exp(-x[i]/epsilon)*(U[i+1] -2*U[i] + U[i-1])/(Math.exp(-x[i+1]/epsilon) -
////                    2*Math.exp(-x[i]/epsilon) + Math.exp(-x[i-1]/epsilon)); //false
//            secondProizvPogranSloiU[i] = (1./(epsilon*epsilon))*Math.exp(-x[i]/epsilon)*(U[i+1] - 2*U[i] + U[i-1])/(Math.exp(-x[i+1]/epsilon) -
//                    2*Math.exp(-x[i]/epsilon) + Math.exp(-x[i-1]/epsilon));
////            System.out.println(proizvPogranSloiU[i]);
//        }

/*

    //вычисление второй производной в 3 точках:
        for (int i=1;i<=L-1;i=i+2){
            secondProizvU[i-1] = (U[i+1] -2*U[i] + U[i-1])/(hh*hh);
            secondProizvU[i] = (U[i+1] -2*U[i] + U[i-1])/(hh*hh);
            secondProizvU[i+1] = (U[i+1] -2*U[i] + U[i-1])/(hh*hh);
        }
        for (int i=1;i<=L-1;i=i+2){
            secondProizvPogranSloiU[i-1] = (U[i-1] - 2.*U[i] + U[i+1]) * (Math.exp(-x[i-1]/epsilon)/(epsilon*epsilon))/(Fi[i-1]-2.*Fi[i] + Fi[i+1]);
            secondProizvPogranSloiU[i] = (U[i-1] - 2.*U[i] + U[i+1]) * (Math.exp(-x[i]/epsilon)/(epsilon*epsilon))/(Fi[i-1]-2.*Fi[i] + Fi[i+1]);
            secondProizvPogranSloiU[i+1] = (U[i-1] - 2.*U[i] + U[i+1]) * (Math.exp(-x[i+1]/epsilon)/(epsilon*epsilon))/(Fi[i-1]-2.*Fi[i] + Fi[i+1]);
        }
*/

        //вычисление второй производной на 5 точках
        for (int i=2;i<=L-2;i=i+4){
            secondProizvU[i-2] = (U[i+1] -2*U[i] + U[i-1])/(hh*hh);
            secondProizvU[i-1] = (U[i+1] -2*U[i] + U[i-1])/(hh*hh);
            secondProizvU[i] = (U[i+1] -2*U[i] + U[i-1])/(hh*hh);
            secondProizvU[i+1] = (U[i+1] -2*U[i] + U[i-1])/(hh*hh);
            secondProizvU[i+2] = (U[i+1] -2*U[i] + U[i-1])/(hh*hh);
        }
        for (int i=2;i<=L-2;i=i+4){
            secondProizvPogranSloiU[i-2] = (U[i-1] - 2.*U[i] + U[i+1]) * (Math.exp(-x[i-2]/epsilon)/(epsilon*epsilon))/(Fi[i-1]-2.*Fi[i] + Fi[i+1]);
            secondProizvPogranSloiU[i-1] = (U[i-1] - 2.*U[i] + U[i+1]) * (Math.exp(-x[i-1]/epsilon)/(epsilon*epsilon))/(Fi[i-1]-2.*Fi[i] + Fi[i+1]);
            secondProizvPogranSloiU[i] = (U[i-1] - 2.*U[i] + U[i+1]) * (Math.exp(-x[i]/epsilon)/(epsilon*epsilon))/(Fi[i-1]-2.*Fi[i] + Fi[i+1]);
            secondProizvPogranSloiU[i+1] = (U[i-1] - 2.*U[i] + U[i+1]) * (Math.exp(-x[i+1]/epsilon)/(epsilon*epsilon))/(Fi[i-1]-2.*Fi[i] + Fi[i+1]);
            secondProizvPogranSloiU[i+2] = (U[i-1] - 2.*U[i] + U[i+1]) * (Math.exp(-x[i+2]/epsilon)/(epsilon*epsilon))/(Fi[i-1]-2.*Fi[i] + Fi[i+1]);
        }

        double maxNorm = 0;
        double[] norm = new double[L+1];
        for (int i=0;i<=L;i++){
//            norm[i] = epsilon*epsilon*Math.abs(secondProizvPogranSloiU[i] - (-Math.PI*Math.PI*Math.cos(Math.PI*x[i]) +
//                    + Math.exp(-x[i]/epsilon)/(epsilon*epsilon)));
            norm[i] = epsilon*epsilon*Math.abs(secondProizvU[i] - (-Math.PI*Math.PI*Math.cos(Math.PI*x[i]) +
                    + Math.exp(-x[i]/epsilon)/(epsilon*epsilon)));

        }

        for(int i=0;i<=L;i++){
            if (maxNorm<norm[i]) maxNorm = norm[i];
        }
        return maxNorm;
    }

    public static double findSecondUsqrt (int N, double epsilon) {
        int L = 3*N;
//        System.out.println(L);
        double[] U = new double[L+1];
        double[] Fi = new double[L+1];
        double[] secondProizvU = new double[L+1];
        double[] secondProizvPogranSloiU = new double[L+1];
//        double h = 1/N;
        double hh = 1./L;
        double[] uzelX = new double[N+1];
        double[] x = new double[L+1];

        x[0] = 0;
        for (int j=1;j<=L;j++){
            x[j] = x[j-1]+hh;
        }

        for (int i=0;i<=L;i++){
            U[i] = Math.cos(Math.PI*x[i]) + Math.exp(-x[i]/epsilon);
//            Fi[i] = Math.exp(-x[i]/epsilon);
            Fi[i] = Math.sqrt(x[i] + epsilon);
        }
        for (int i=1;i<=L-1;i=i+3){
            secondProizvU[i-1] = (U[i+1] -2*U[i] + U[i-1])/(hh*hh);
            secondProizvU[i] = (U[i+1] -2*U[i] + U[i-1])/(hh*hh);
            secondProizvU[i+1] = (U[i+1] -2*U[i] + U[i-1])/(hh*hh);
        }
        for (int i=1;i<=L-1;i=i+3){
//            secondProizvPogranSloiU[i] = (1./(epsilon*epsilon)) * Math.exp(-x[i]/epsilon)*(U[i+1] - 2*U[i] + U[i-1])/(Math.exp(-x[i+1]/epsilon) -
//                    2*Math.exp(-x[i]/epsilon) + Math.exp(-x[i-1]/epsilon));
            secondProizvPogranSloiU[i-1] = (U[i-1] - 2.*U[i] + U[i+1]) * (-1./(4.*(Math.sqrt((x[i-1]+epsilon)*(x[i-1]+epsilon)*(x[i-1]+epsilon)))))/(Fi[i-1]-2.*Fi[i] + Fi[i+1]);
            secondProizvPogranSloiU[i] = (U[i-1] - 2.*U[i] + U[i+1]) * (-1./(4.*(Math.sqrt((x[i]+epsilon)*(x[i]+epsilon)*(x[i]+epsilon)))))/(Fi[i-1]-2.*Fi[i] + Fi[i+1]);
            secondProizvPogranSloiU[i+1] = (U[i-1] - 2.*U[i] + U[i+1]) * (-1./(4.*(Math.sqrt((x[i+1]+epsilon)*(x[i+1]+epsilon)*(x[i+1]+epsilon)))))/(Fi[i-1]-2.*Fi[i] + Fi[i+1]);
        }
        double maxNorm = 0;
        double[] norm = new double[L+1];
        for (int i=0;i<=L;i++){
            norm[i] = epsilon*epsilon*Math.abs(secondProizvPogranSloiU[i] - (-Math.PI*Math.PI*Math.cos(Math.PI*x[i]) +
                    + Math.exp(-x[i]/epsilon)/(epsilon*epsilon)));
//            norm[i] = epsilon*epsilon*Math.abs(secondProizvU[i] - (-Math.PI*Math.PI*Math.cos(Math.PI*x[i]) +
//                    + Math.exp(-x[i]/epsilon)/(epsilon*epsilon)));

        }

        for(int i=0;i<=L-1;i++){
            if (maxNorm<norm[i]) maxNorm = norm[i];
        }
        return maxNorm;
    }



    public static void main(String[] args) {
        double a = 0.;
        double epsilon = 1./512.;
        for (int i=16;i<=1024;i=2*i){
//            System.out.println(findU(i, 1.));
            double b = findUexp(i, epsilon);
            double four = a/b;
            a = b;
            System.out.println("epsilon = 1/"+1./epsilon);
            System.out.println("Порядок точности log2 (||"+i/2.+"||/||"+i+"||) = "+Math.log10(four)/Math.log10(2.));
            System.out.println("_______________________________________________________________");
            System.out.println("||"+i+"|| = "+ b);
        }
    }
}
