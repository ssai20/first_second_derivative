import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;

enum CodeError
{
    WRONG_CONSOLE_CONDITIONS("Введите, пожалуйста, данные по порядку: (вид сортировки) (тип сортируемых данных) (путь к исходящему файлу в формате txt) (пути ко входящим файлам в форматах txt)."),
    EMPTY_FIELD_TYPE("Введите, пожалуйста, тип сортируемых данных в формате '-i' или '-s'!"),
    FILE_ERROR("Извините, файла %s не существует, либо указан неверный путь."),
    EMPTY_FILE("Введите пути к файлам!"),
    WRONG_TYPE_INTEGER("Вы ввели '-i' в консоль, но %s строка в файле %s не типа Integer, либо число выходит за рамки диапазона Integer"), //Минимальное значение числа типа int: -2 147 483 648, максимальное значение: 2 147 483 647
    WRONG_SPACE("В строках файла не могут быть пробелы, но %s строка в файле %s имеет пробел");

    private String errorString;
    CodeError(String error) {
        this.errorString = error;
    }

    public String getErrorString() {
        return errorString;
    }

}


class CodeException extends Exception
{
    private CodeError codeError;
    public CodeException (CodeError codeError)
    {
        super(codeError.getErrorString());
    }
    public CodeException (CodeError codeError, String param)
    {
        super (String.format(codeError.getErrorString(),param));
    }

    public CodeException (CodeError codeError, String param1, String param2)
    {
        super (String.format(codeError.getErrorString(),param1, param2));
    }

    public CodeError getCodeError() {
        return codeError;
    }
}

class GenericTypes<T extends Comparable<T>> implements Comparable<T>{
    public T a;

    public boolean isListConsistence(List<T> list){
            boolean isConsistence = true;
            if (list.size()==1){
                return true;
            }
            for (int i=0;i<list.size()-1;i++){
                if (list.get(i).compareTo(list.get(i+1))>0){
                    return false;
                }
            }
            return isConsistence;
    }

    public boolean isNumeric(String str) {
        try {
            Integer.parseInt(str);
            return true;
        } catch (NumberFormatException e) {
            return false;
        }
    }
    public void successFinish(List<T> list){
        if (isListConsistence(list)) {
            System.out.println("Сортировка успешно завершена!");
            }
    }


    public List<T> genericMergeSortAlgorithm(List<T> listFile1, List<T> listFile2) {

        if (!isListConsistence(listFile1)){
            List<T> deepList1 = new ArrayList<T>();
            List<T> deepList2 = new ArrayList<T>();
            for (int i=0;i<(listFile1.size()/2);i++){
                deepList1.add(listFile1.get(i));
            }
            for (int i=listFile1.size()/2;i<listFile1.size();i++){
                deepList2.add(listFile1.get(i));
            }
            listFile1 = genericMergeSortAlgorithm(deepList1, deepList2);
        }

        if (!isListConsistence(listFile2)){
            List<T> deepList3 = new ArrayList<T>();
            List<T> deepList4 = new ArrayList<T>();
            for (int i=0;i<(listFile2.size()/2);i++){
                deepList3.add(listFile2.get(i));
            }
            for (int i=listFile2.size()/2;i<listFile2.size();i++){
                deepList4.add(listFile2.get(i));
            }
            listFile2 = genericMergeSortAlgorithm(deepList3, deepList4);
        }

        List<T> listResult = new ArrayList<T>();
        if ((Objects.equals(listFile1, new ArrayList<T>()))&&(!Objects.equals(listFile2, new ArrayList<T>()))){
            return listFile2;
        }
        if ((!Objects.equals(listFile1, new ArrayList<T>()))&&(Objects.equals(listFile2, new ArrayList<T>()))){
            return listFile1;
        }


        Queue<T> queueFile1 = new LinkedList<T>(listFile1);
        Queue<T> queueFile2 = new LinkedList<T>(listFile2);
        T buffer = null;
        boolean bufferInRightColumn = true;
        if (queueFile1.peek().compareTo(queueFile2.peek()) < 0) {
            buffer = queueFile2.peek();
            bufferInRightColumn = true;
            listResult.add(queueFile1.peek());
            queueFile1.poll();
        } else if (queueFile1.peek().compareTo(queueFile2.peek()) >= 0) {
            bufferInRightColumn = false;
            buffer = queueFile1.peek();
            listResult.add(queueFile2.peek());
            queueFile2.poll();
        }

        while (!queueFile1.isEmpty() || !queueFile2.isEmpty()) {
            if (!queueFile1.isEmpty() && !queueFile2.isEmpty()) {
                if (bufferInRightColumn) {
                    if (queueFile1.peek().compareTo(buffer) < 0) {
                        listResult.add(queueFile1.peek());
                        queueFile1.poll();
//                        bufferInRightColumn = true;
                    } else if (queueFile1.peek().compareTo(buffer) >= 0) {
                        listResult.add(buffer);
                        buffer = queueFile1.peek();
                        queueFile2.poll();
                        bufferInRightColumn = false;
                    }

                } else {
                    if (queueFile2.peek().compareTo(buffer) < 0) {
                        listResult.add(queueFile2.peek());
                        queueFile2.poll();
                        bufferInRightColumn = false;
                    } else if (queueFile2.peek().compareTo(buffer) >= 0){
                        listResult.add(buffer);
                        buffer = queueFile2.peek();
                        queueFile1.poll();
                        bufferInRightColumn = true;
                    }
                }
            }


            if (queueFile1.isEmpty() && !queueFile2.isEmpty()) {
                listResult.add(buffer);
                queueFile2.poll();
                listResult.addAll(queueFile2);
                queueFile2.removeAll(queueFile2);
            }
            if (!queueFile1.isEmpty() && queueFile2.isEmpty()) {
                listResult.add(buffer);
                queueFile1.poll();
                listResult.addAll(queueFile1);
                queueFile1.removeAll(queueFile1);
            }
        }
        return listResult;

    }



    public List<T> mapFileToList(String fileInPath, String dataType) throws CodeException {
        if (!Files.exists(Path.of(fileInPath))) throw new CodeException(CodeError.FILE_ERROR, fileInPath);
        File file = new File(fileInPath);
        List<T> listFile = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(file), "UTF-8"))) {

                int stringNumber = 0;
                while (true) {
                    stringNumber++;
                    String string = br.readLine();
                    if (string == null) {
                        break;
                    }


                    if (string.contains(" "))
                        throw new CodeException(CodeError.WRONG_SPACE, String.valueOf(stringNumber), file.toString());
                    if (dataType.equals("-i")) {
                        if (!isNumeric(string))
                            throw new CodeException(CodeError.WRONG_TYPE_INTEGER, String.valueOf(stringNumber), file.toString());
                        Integer i = Integer.parseInt(string);
                        listFile.add((T) i);
                    }

                    if (dataType.equals("-s")) {
                        listFile.add((T) string);
                    }
                }
            } catch (IOException e1) {
                e1.printStackTrace();
            }

            if (!isListConsistence(listFile)) {
                System.out.println("Файл " + fileInPath + " не отcортирован по возрастанию, сортировка будет продолжена");
            }

        return listFile;
    }

    public void mapListToFile (String fileOutPath, List<T> resultList) throws CodeException {
        File file = new File(fileOutPath);
        try (BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(file),  "UTF-8"))){
            for (int i=0;i<resultList.size();i++) {
                bw.write(resultList.get(i).toString());
                bw.newLine();
            }
        } catch (IOException e1){
            e1.printStackTrace();
        }
        System.out.println("Данные отправлены в файл: "+fileOutPath);
    }

    @Override
    public int compareTo(T o) {
        if ((this.a instanceof Integer) && (o instanceof Integer)){
            return (Integer)this.a-(Integer)o;
        }
        if ((this.a instanceof String) && (o instanceof String)){
            return ((String) this.a).compareTo((String) o);
        }
        return 0;
    }
}




public class SortIt{
    public static void mergeSort(String conditions) throws CodeException {
        try {

        String[] s = conditions.split(" ");
        if (!s[0].equals("-a") && !s[0].equals("-d")) {
            StringBuilder conditionsSB = new StringBuilder(conditions);
            conditionsSB.insert(0, "-a ");
            conditions = conditionsSB.toString();
        }
        List<String> conditionList = List.of(conditions.split(" "));
        if ((!conditionList.get(1).equals("-s")) && (!conditionList.get(1).equals("-i"))) throw new CodeException(CodeError.EMPTY_FIELD_TYPE);
        if ((conditionList.size()==2)) throw new CodeException(CodeError.EMPTY_FIELD_TYPE);
        if (!Files.exists(Path.of(conditionList.get(2)))) throw new CodeException(CodeError.FILE_ERROR, conditionList.get(2));

        String sortType = conditionList.get(0);
        String dataType = conditionList.get(1);
        String fileOutPath = conditionList.get(2);


        if (dataType.equals("-i")) {
            GenericTypes<Integer> genericTypes = new GenericTypes<>();
            List<Integer> resultList = new ArrayList<>();

            for (int i = 3; i < conditionList.size(); i++) {
                resultList = genericTypes.genericMergeSortAlgorithm(resultList, genericTypes.mapFileToList(conditionList.get(i), dataType));
            }
            if (sortType.equals("-a")) {
                genericTypes.mapListToFile(fileOutPath, resultList);
                genericTypes.successFinish(resultList);
            }
            if (sortType.equals("-d")) {
                genericTypes.successFinish(resultList);
                Collections.reverse(resultList);
                genericTypes.mapListToFile(fileOutPath, resultList);
            }
        }

        if (dataType.equals("-s")) {
            GenericTypes<String> genericTypes = new GenericTypes<>();
            List<String> resultList = new ArrayList<>();
            for (int i = 3; i < conditionList.size(); i++) {
                resultList = genericTypes.genericMergeSortAlgorithm(resultList, genericTypes.mapFileToList(conditionList.get(i), dataType));
            }
            if (sortType.equals("-a")) {
                genericTypes.mapListToFile(fileOutPath, resultList);
                genericTypes.successFinish(resultList);
            }
            if (sortType.equals("-d")) {
                genericTypes.successFinish(resultList);
                Collections.reverse(resultList);
                genericTypes.mapListToFile(fileOutPath, resultList);
            }
        }
    } catch (CodeException c){
        System.out.println(c.getMessage());
    }

    }
    public static void main(String[] args) throws CodeException {
        try {
            String conditions = null;
            try (BufferedReader br = new BufferedReader(new InputStreamReader(System.in))) {
                conditions = br.readLine();
            } catch (IOException e) {
                e.printStackTrace();
            }

            if ((conditions.equals("-i"))||(conditions.equals("-s"))) throw new CodeException(CodeError.EMPTY_FILE);
            if (conditions.equals("")||!conditions.contains(" ")||conditions.matches("^\\s*$")) throw new CodeException(CodeError.WRONG_CONSOLE_CONDITIONS);

            mergeSort(conditions) ;
        } catch (CodeException c) {
            System.out.println(c.getMessage());
        }
    }
}


//-d -i /Users/work/Desktop/out.txt /Users/work/Desktop/in1.txt /Users/work/Desktop/in2.txt



