import java.util.*;

public class Sculptures {

    public static void idealSculpturesAndTheirIndexes(int N, long X, long T, int[] sculptureWeights)
    {

        if (N < 1 || N > 200000) {
            throw new RuntimeException("Недопустимое значение для N");
        }

        if (X < 0 || X > 1000000000) {
            throw new RuntimeException("Недопустимое значение для X");
        }

        if (T < 0 || T > 300000000000000L) {
            throw new RuntimeException("Недопустимое значение для T");
        }



        TreeMap<Long, Integer> mapMinuteToIndex = new TreeMap();
        for (int i=0;i<N;i++)
        {
            long weightsDifference = Math.abs(sculptureWeights[i]-X);
            if ((weightsDifference >= 0) && (weightsDifference <= T))
            {
                mapMinuteToIndex.put(weightsDifference, i);
            }
        }


        List<Long> minuteList = new ArrayList<>(mapMinuteToIndex.keySet());
        if (minuteList.size()!=0) {
            long commonWeight = 0;
            int k = 0;
            for (int i=0;i<minuteList.size();i++)
            {
                commonWeight = commonWeight + minuteList.get(i);
                if (commonWeight<=T) {
                    k++;
                }
            }


            System.out.println(k);
            int a = mapMinuteToIndex.get(minuteList.get(0))+1;
            System.out.print(a);
            if (k>1) {
                for (int i = 1; i < k; i++) {
                    int b = mapMinuteToIndex.get(minuteList.get(i))+1;
                    System.out.print(" " + b);
                }
            }
        }

        if (minuteList.size()==0)
        {
            System.out.println(0);
        }
    }

    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);

        int N = scanner.nextInt();

        long X = scanner.nextLong();

        long T = scanner.nextLong();

        int[] sculptureWeights = new int[N];
        for (int i = 0; i < N; i++) {
            sculptureWeights[i] = scanner.nextInt();
        }

        idealSculpturesAndTheirIndexes(N, X, T, sculptureWeights);
    }
}
