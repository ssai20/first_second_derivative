import java.util.*;

public class ChessGame {

    public static class BlackKing {
        int kingPositionI,kingPositionJ, N;
        public BlackKing(int i, int j, int N)
        {
            this.kingPositionI = i;
            this.kingPositionJ = j;
            this.N = N;
        }

        public int getKingPositionI() {
            return kingPositionI;
        }

        public int getKingPositionJ() {
            return kingPositionJ;
        }

        public void setPosition(int positionI, int positionJ) {
            this.kingPositionI = positionI;
            this.kingPositionJ = positionJ;
        }
    }

    public static class WhitePawn {
        int pawnPositionI, pawnPositionJ, N;
        boolean isExist = false;
        public WhitePawn(int i, int j, int N)
        {
            this.pawnPositionI = i;
            this.pawnPositionJ = j;
            this.N = N;
        }


        public void setExist(boolean exist) {
            isExist = exist;
        }
        public boolean getExist() {
            return isExist;
        }

        public int getPawnPositionI() {
            return pawnPositionI;
        }

        public int getPawnPositionJ() {
            return pawnPositionJ;
        }

        public void setPosition(int positionI, int positionJ) {
            this.pawnPositionI = positionI;
            this.pawnPositionJ = positionJ;
        }
    }

    public static class WhiteRook {
        int rookPositionI, rookPositionJ, N;
        boolean isExist = false;


        public WhiteRook(int i, int j, int N) {
            this.rookPositionI = i;
            this.rookPositionJ = j;
            this.N= N;
        }


        public void setExist(boolean exist) {
            isExist = exist;
        }
        public boolean getExist() {
            return isExist;
        }

        public int getRookPositionI() {
            return rookPositionI;
        }

        public int getRookPositionJ() {
            return rookPositionJ;
        }


        public void setPosition(int positionI, int positionJ) {
            this.rookPositionI = positionI;
            this.rookPositionJ = positionJ;
        }
    }

    public static class WhiteKnight {
        int knightPositionI, knightPositionJ, N;
        boolean isExist = false;
        public WhiteKnight(int i, int j, int N) {
            this.knightPositionI = i;
            this.knightPositionJ = j;
            this.N = N;
        }


        public int getKnightPositionI() {
            return knightPositionI;
        }

        public int getKnightPositionJ() {
            return knightPositionJ;
        }

        public void setExist(boolean exist) {
            isExist = exist;
        }
        public boolean getExist() {
            return isExist;
        }

        public void setPosition(int positionI, int positionJ) {
            this.knightPositionI = positionI;
            this.knightPositionJ = positionJ;
        }
    }

    public static class WhiteBishop {
        int bishopPositionI, bishopPositionJ, N;
        boolean isExist = false;
        public WhiteBishop(int i, int j, int N) {
            this.bishopPositionI = i;
            this.bishopPositionJ = j;
            this.N = N;
        }

        public void setExist(boolean exist) {
            isExist = exist;
        }
        public boolean getExist() {
            return isExist;
        }
        public int getBishopPositionI() {
            return bishopPositionI;
        }

        public int getBishopPositionJ() {
            return bishopPositionJ;
        }

        public void setPosition(int positionI, int positionJ) {
            this.bishopPositionI = positionI;
            this.bishopPositionJ = positionJ;
        }
    }

    public static class WhiteQueen {
        int queenPositionI, queenPositionJ, N;
        boolean isExist = false;

        public WhiteQueen(int i, int j, int N) {
            this.queenPositionI = i;
            this.queenPositionJ = j;
            this.N = N;
        }

        public int getQueenPositionI() {
            return queenPositionI;
        }

        public int getQueenPositionJ() {
            return queenPositionJ;
        }

        public void setExist(boolean exist) {
            isExist = exist;
        }
        public boolean getExist() {
            return isExist;
        }

        public void setPosition(int positionI, int positionJ) {
            this.queenPositionI = positionI;
            this.queenPositionJ = positionJ;
        }
    }


    private static void mapSymbolsToDesk(String area, int N, BlackKing blackKing, WhiteQueen[] whiteQueens, WhiteBishop[] whiteBishops, WhiteKnight[] whiteKnights, WhiteRook[] whiteRooks, WhitePawn[] whitePawns, int[][] chessDesk)
    {

        int k = 0;
        int qCounter = 0;
        int bCounter = 0;
        int knCounter = 0;
        int rCounter = 0;
        int pCounter = 0;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (area.charAt(k) == 'K')
                {
                    blackKing.setPosition(i,j);
                }
                if (area.charAt(k) == 'Q')
                {
                    whiteQueens[qCounter].setPosition(i,j);
                    whiteQueens[qCounter].setExist(true);
                    qCounter++;
                }
                if (area.charAt(k) == 'B')
                {
                    whiteBishops[bCounter].setPosition(i,j);
                    whiteBishops[bCounter].setExist(true);
                    bCounter++;
                }
                if (area.charAt(k) == 'N')
                {
                    whiteKnights[knCounter].setPosition(i,j);
                    whiteKnights[knCounter].setExist(true);
                    knCounter++;
                }
                if (area.charAt(k) == 'R')
                {
                    whiteRooks[rCounter].setPosition(i,j);
                    whiteRooks[rCounter].setExist(true);
                    rCounter++;
                }
                if (area.charAt(k) == 'P')
                {
                    whitePawns[pCounter].setPosition(i,j);
                    whitePawns[pCounter].setExist(true);
                }

                k++;
            }
        }
        for (int i=0;i<N;i++)
        {
            for (int j = 0;j<N;j++)
            {
                chessDesk[i][j] = 1;
            }
        }
        for (WhiteQueen whiteQueen:whiteQueens) {
            if (whiteQueen.getExist()) {
                chessDesk[whiteQueen.getQueenPositionI()][whiteQueen.getQueenPositionJ()] = 50;
            }
        }
        for (WhiteBishop whiteBishop:whiteBishops) {
            if (whiteBishop.getExist()) {
                chessDesk[whiteBishop.getBishopPositionI()][whiteBishop.getBishopPositionJ()] = 51;
            }
        }
        for (WhiteRook whiteRook:whiteRooks) {
            if (whiteRook.getExist()) {
                chessDesk[whiteRook.getRookPositionI()][whiteRook.getRookPositionJ()] = 52;
            }
        }
        for (WhiteKnight whiteKnight:whiteKnights) {
            if (whiteKnight.getExist()) {
                chessDesk[whiteKnight.getKnightPositionI()][whiteKnight.getKnightPositionJ()] = 53;
            }
        }
        for (WhitePawn whitePawn:whitePawns) {
            if (whitePawn.getExist()) {
                chessDesk[whitePawn.getPawnPositionI()][whitePawn.getPawnPositionJ()] = 50;
            }
        }

    }

    public static boolean isBlackKingUnderAttack(BlackKing blackKing, WhiteQueen[] whiteQueens, WhiteBishop[] whiteBishops, WhiteKnight[] whiteKnights, WhiteRook[] whiteRooks, WhitePawn[] whitePawns, int[][] chessDesk)
    {
        boolean isKnK = false;
        boolean isPK = false;
        boolean isQK = false;
        boolean isBK = false;
        boolean isRK = false;
        for (WhiteKnight whiteKnight:whiteKnights) {
            isKnK = isKnK||isKnightAttackKing(blackKing.getKingPositionI(), blackKing.getKingPositionJ(), whiteKnight.getKnightPositionI(), whiteKnight.getKnightPositionJ()) && whiteKnight.getExist();
        }
        for (WhitePawn whitePawn:whitePawns) {
            isPK = isPK || isPawnAttackKing(blackKing.getKingPositionI(), blackKing.getKingPositionJ(), whitePawn.getPawnPositionI(), whitePawn.getPawnPositionJ()) && whitePawn.getExist();
        }

        for (WhiteQueen whiteQueen:whiteQueens) {
            isQK = isQK||isQueenAttackKing(blackKing.getKingPositionI(), blackKing.getKingPositionJ(), whiteQueen.getQueenPositionI(), whiteQueen.getQueenPositionJ(), chessDesk) && whiteQueen.getExist();
        }
        for (WhiteBishop whiteBishop:whiteBishops) {
            isBK = isBK||isBishopAttackKing(blackKing.getKingPositionI(), blackKing.getKingPositionJ(), whiteBishop.getBishopPositionI(), whiteBishop.getBishopPositionJ(), chessDesk) && whiteBishop.getExist();
        }
        for (WhiteRook whiteRook:whiteRooks) {
            isRK = isRK||isRookAttackKing(blackKing.getKingPositionI(), blackKing.getKingPositionJ(), whiteRook.getRookPositionI(), whiteRook.getRookPositionJ(), chessDesk) && whiteRook.getExist();
        }

        System.out.println("Queen = "+isQK);
        System.out.println("Bishop = "+isBK);
        System.out.println("Knight = "+isKnK);
        System.out.println("Rook = "+isRK);
        System.out.println("Pawn = "+isPK);
        return isKnK||isPK||isQK||isBK||isRK;
    }

    public static int lengthBetweenQueenKing(int iK, int jK, int iQ, int jQ)
    {
        int difference = 190;

        if (iK==iQ) difference = Math.abs(jK-jQ);
        if (jK==jQ) difference = Math.abs(iK-iQ);
        if (Math.abs(iK - iQ) == Math.abs(jK - jQ) && jK != jQ) difference = Math.abs(iK-iQ);
        return difference;
    }

    public static int lengthBetweenRookKing(int iK, int jK, int iR, int jR)
    {
        int difference = 180;

        if (iK==iR) difference = Math.abs(jK-jR);
        if (jK==jR) difference = Math.abs(iK-iR);
        return difference;
    }
    public static int lengthBetweenBishopKing(int iK, int jK, int iB, int jB)
    {
        int difference = 170;

        if (Math.abs(iK-iB)==Math.abs(jK-jB)&& jK != jB)
        {
            difference = Math.abs(iK-iB);
        }
        return difference;
    }
    public static boolean isPawnAttackKing(int iK, int jK, int iP, int jP)
    {
        return (iP - iK == 1) && (Math.abs(jP - jK) == 1);
    }

    public static boolean isKnightAttackKing(int iK, int jK, int iP, int jP)
    {
        return ((Math.abs(iK - iP) == 2) && (Math.abs(jK - jP) == 1)) || ((Math.abs(jK - jP) == 2) && (Math.abs(iK - iP) == 1));
    }

    public static boolean isQueenAttackKing(int iK, int jK, int iQ, int jQ, int[][] matrix)
    {
        boolean isQueenAttack = false;
        int difference = 230;
        if ((iK==iQ)&&(jK<jQ))
        {
            difference = 0;
            for (int i=jK;i<jQ;i++)
            {
                difference = difference + matrix[iK][i];
            }
        }

        if ((iK==iQ)&&(jK>jQ))
        {
            difference = 0;
            for (int i=jK;i>jQ;i--)
            {
                difference = difference + matrix[iK][i];
            }
        }
        if ((jK==jQ)&&(iK>iQ))
        {
            difference = 0;
            for (int i=iK;i>iQ;i--)
            {
                difference = difference + matrix[i][jK];
            }
        }

        if ((jK==jQ)&&(iK<iQ))
        {
            difference = 0;
            for (int i=iK;i<iQ;i++)
            {
                difference = difference + matrix[i][jK];
            }
        }
        if (Math.abs(iK-iQ)==Math.abs(jK-jQ)&&(iK>iQ)&&(jK>jQ))
        {
            difference = 0;
            for (int i=iK,j=jK;i>iQ;i--)
            {
                difference = difference + matrix[i][j];
                j--;
            }
        }

        if ((Math.abs(iK-iQ)==Math.abs(jK-jQ))&&(iK<iQ)&&(jK<jQ))
        {
            difference = 0;
            for (int i=iK, j=jK;i<iQ;i++)
            {
                difference = difference + matrix[i][j];
                j++;
            }
        }
        if ((Math.abs(iK-iQ)==Math.abs(jK-jQ))&&(iK<iQ)&&(jK>jQ))
        {
            difference = 0;
            for (int i=iK,j=jK;i<iQ;i++)
            {
                difference = difference + matrix[i][j];
                j--;
            }
        }

        if ((Math.abs(iK-iQ)==Math.abs(jK-jQ))&&(iK>iQ)&&(jK<jQ))
        {
            difference = 0;
            for (int i=iK,j=jK;i>iQ;i--)
            {
                difference = difference + matrix[i][j];
                j++;
            }
        }
        if (difference==lengthBetweenQueenKing(iK,jK,iQ,jQ)) isQueenAttack = true;
        return isQueenAttack;
    }

    public static boolean isBishopAttackKing(int iK, int jK, int iQ, int jQ, int[][] matrix)
    {
        boolean isBishopAttack = false;
        int difference = 120;

        if ((Math.abs(iK-iQ)==Math.abs(jK-jQ))&&(iK>iQ)&&(jK>jQ))
        {
            difference = 0;
            for (int i=iK,j=jK;i>iQ;i--)
            {
                difference = difference + matrix[i][j];
                j--;
            }
        }

        if ((Math.abs(iK-iQ)==Math.abs(jK-jQ))&&(iK<iQ)&&(jK<jQ))
        {
            difference = 0;
            for (int i=iK, j=jK;i<iQ;i++)
            {
                difference = difference + matrix[i][j];
                j++;
            }
        }
        if ((Math.abs(iK-iQ)==Math.abs(jK-jQ))&&(iK<iQ)&&(jK>jQ))
        {
            difference = 0;
            for (int i=iK,j=jK;i<iQ;i++)
            {
                difference = difference + matrix[i][j];
                j--;
            }
        }

        if ((Math.abs(iK-iQ)==Math.abs(jK-jQ))&&(iK>iQ)&&(jK<jQ))
        {
            difference = 0;
            for (int i=iK,j=jK;i>iQ;i--)
            {
                difference = difference + matrix[i][j];
                j++;
            }
        }
        if (difference==lengthBetweenBishopKing(iK,jK,iQ,jQ)) isBishopAttack = true;
        return isBishopAttack;
    }

    public static boolean isRookAttackKing(int iK, int jK, int iQ, int jQ, int[][] matrix)
    {
        boolean isRookAttack = false;
        int difference = 1000;
        if ((iK==iQ)&&(jK<jQ))
        {
            difference = 0;
            for (int i=jK;i<jQ;i++)
            {
                difference = difference + matrix[iK][i];
            }
        }

        if ((iK==iQ)&&(jK>jQ))
        {
            difference = 0;
            for (int i=jK;i>jQ;i--)
            {
                difference = difference + matrix[iK][i];
            }
        }

        if ((jK==jQ)&&(iK>iQ))
        {
            difference = 0;
            for (int i=iK;i>iQ;i--)
            {
                difference = difference + matrix[i][jK];
            }
        }

        if ((jK==jQ)&&(iK<iQ))
        {
            difference = 0;
            for (int i=iK;i<iQ;i++)
            {
                difference = difference + matrix[i][jK];
            }
        }

        if ((difference==lengthBetweenRookKing(iK,jK,iQ,jQ))) isRookAttack = true;
        return isRookAttack;
    }

    public static boolean isKingAttacked(String area, int deskRank)
    {
        if ( (area.length()>16) || (area.length()<16) ) throw new RuntimeException("Path should have length equals 16");
        boolean isKingAttacked;
        int[][] chessDesk = new int[deskRank][deskRank];
        BlackKing blackKing = new BlackKing(0,0, deskRank);

        WhiteQueen[] whiteQueens = new WhiteQueen[deskRank*deskRank];
        for (int i = 0; i < whiteQueens.length; i++) {
            whiteQueens[i] = new WhiteQueen(0, 0, deskRank);
        }

        WhiteBishop[] whiteBishops = new WhiteBishop[deskRank*deskRank];
        for (int i = 0; i < whiteBishops.length; i++) {
            whiteBishops[i] = new WhiteBishop(0, 0, deskRank);
        }

        WhiteKnight[] whiteKnights = new WhiteKnight[deskRank*deskRank];
        for (int i = 0; i < whiteKnights.length; i++) {
            whiteKnights[i] = new WhiteKnight(0, 0, deskRank);
        }

        WhiteRook[] whiteRooks = new WhiteRook[deskRank*deskRank];
        for (int i = 0; i < whiteRooks.length; i++) {
            whiteRooks[i] = new WhiteRook(0, 0, deskRank);
        }

        WhitePawn[] whitePawns = new WhitePawn[deskRank*deskRank];
        for (int i = 0; i < whitePawns.length; i++) {
            whitePawns[i] = new WhitePawn(0, 0, deskRank);
        }

        mapSymbolsToDesk(area, deskRank, blackKing, whiteQueens, whiteBishops, whiteKnights, whiteRooks, whitePawns, chessDesk);
        isKingAttacked = isBlackKingUnderAttack(blackKing, whiteQueens, whiteBishops, whiteKnights, whiteRooks, whitePawns, chessDesk);
        return isKingAttacked;
    }

    public static void main(String[] args) {
        Scanner scan = new Scanner(System.in);
        String area = scan.nextLine();
        boolean isKingAttacked = isKingAttacked(area, 4);
        System.out.println(isKingAttacked);
    }
}
