import java.util.*;

public class Player {
    private static final int NEGINFINITY = Integer.MIN_VALUE;
    private static final int POSINFINITY = Integer.MAX_VALUE;
    private final long LIMIT =  410000000L; //timelmit for each board is max of 0.5 seconds. 
    private long TIME;
    private final int [][][] MDIAGS = { {{0,0,0},{1,1,1},{2,2,2},{3,3,3}}, {{3,0,0},{2,1,1},{1,2,2},{0,3,3}}, 
                                        {{0,3,0},{1,2,1},{2,1,2},{3,0,3}}, {{3,3,0},{2,2,1},{1,1,2},{0,0,3}}  };

    /**
     * Performs a move
     *
     * @param gameState
     *            the current state of the board
     * @param deadline
     *            time before which we must have returned
     * @return the next state the board is in after our move
     */
  public GameState play(final GameState gameState, final Deadline deadline) {
        Vector<GameState> nextStates = new Vector<GameState>();
        TIME = deadline.timeUntil();
        gameState.findPossibleMoves(nextStates);
        int vSize = nextStates.size();

        if (nextStates.size() == 0) {
            // Must play "pass" move if there are no other moves possible.
            return new GameState(gameState, new Move());
        }

        int argmax=-1;
      
        int bestValue = NEGINFINITY;
        int maxDepth =3;

            int alpha =  Integer.MIN_VALUE;
            int beta =  Integer.MAX_VALUE;
            for (int i =0; i < vSize; i++)
             {
                int res=0;
                if (gameState.getNextPlayer() == Constants.CELL_O){
                   res = minValue(nextStates.get(i),1,maxDepth, alpha, beta, deadline);
                }
                 if (gameState.getNextPlayer() == Constants.CELL_X){
                   res = maxValue(nextStates.get(i),1,maxDepth, alpha, beta, deadline);
                }
                
                if (res >= bestValue){
                    bestValue = res;
                    argmax = i;
                }
  
        } 
    
    return nextStates.get(argmax);     

    }
    private int minValue(GameState gState, int depth, int maxDepth, int alpha, int beta, Deadline deadline) {

         if (gState.isEOG()) {

            return evalTerminalState(gState);
         }

         if (depth == maxDepth) {
            return evalState(gState);
         }


        Vector<GameState> nextStates = new Vector<GameState>();
        gState.findPossibleMoves(nextStates);

        int v = POSINFINITY;
        for (GameState state : nextStates){

            v = Math.min(v, maxValue(state, depth+1, maxDepth, alpha, beta, deadline));
            if (v <= alpha){
                return v;
            }
            beta = Math.min(beta,v);
        }
        return v; 
        
    }

    private int maxValue(GameState gState, int depth, int maxDepth, int alpha, int beta, Deadline deadline) {

         if (gState.isEOG()) { // must be a terminal state
    
            return evalTerminalState(gState);
         }

         if (depth == maxDepth) {
            return evalState(gState);
         }
    
        Vector<GameState> nextStates = new Vector<GameState>();
        gState.findPossibleMoves(nextStates);

        int v = NEGINFINITY;
        for (GameState state : nextStates){

            v = Math.max(v, minValue(state, depth+1, maxDepth, alpha, beta, deadline));
            if (v >= beta){
                return v;
            }
            alpha = Math.max(alpha,v);
        }
        return v; 
    }

    private int evalState(final GameState gameState) 
    {

        final int BOARD_SIZE = gameState.BOARD_SIZE;
        int nrRowsWinnableX =0; 
        int nrRowsWinnableO =0;   
        int nrColWinnableX =0; 
        int nrColWinnableO =0; 
        int nrDiagWinnableX =0; 
        int nrDiagWinnableO =0;
        int nrOrthWinnableX =0; 
        int nrOrthWinnableO =0;
        int zDiagX =0; 
        int zDiagO =0;
        int yDiagX =0; 
        int yDiagO =0;
        int xDiagX =0; 
        int xDiagO =0;

        for(int layer =0; layer < BOARD_SIZE; layer++) {

            for (int row=0; row < BOARD_SIZE; row++) {

              int possibleWinX =0; int possibleWin0 =0;

                for (int col=0; col < BOARD_SIZE; col++) {

                   int cell = gameState.at(row,col,layer);

                   if (cell == Constants.CELL_EMPTY){
                        possibleWinX++;
                        possibleWin0++;
                   }

                   if (cell == Constants.CELL_X)           
                        possibleWinX++;

                   if(cell == Constants.CELL_O)
                        possibleWin0++;

                } // end inner for loop

                if (possibleWinX == BOARD_SIZE)
                    nrRowsWinnableX++;
                if (possibleWin0 == BOARD_SIZE)
                    nrRowsWinnableO++;

            } 
         } // end outer for loop

        // Process the columns layer by layer
        for(int layer =0; layer < BOARD_SIZE; layer++) {

            for (int col =0; col< BOARD_SIZE; col++) {

                int possibleWinX =0; int possibleWin0 =0;

                for (int row=0; row < BOARD_SIZE; row++) {

                   int cell = gameState.at(row,col, layer);

                   if (cell == Constants.CELL_EMPTY){
                        possibleWinX++;
                        possibleWin0++;
                   }
                   if (cell == Constants.CELL_X)           
                        possibleWinX++;

                   if(cell == Constants.CELL_O)
                        possibleWin0++;

                } // end inner for loop

                if (possibleWinX == BOARD_SIZE)
                    nrColWinnableX++;
                if (possibleWin0 == BOARD_SIZE)
                    nrColWinnableO++;

            } 
        }// end outer for loop

      // Process the orthogonal wins through all the layers
     for (int row=0; row < BOARD_SIZE; row++) {

          for (int col=0; col < BOARD_SIZE; col++) {

              int possibleWinX =0; int possibleWin0 =0;

                for(int layer =0; layer < BOARD_SIZE; layer++)  {

                   int cell = gameState.at(row,col,layer);

                   if (cell == Constants.CELL_EMPTY){
                        possibleWinX++;
                        possibleWin0++;
                   }

                   if (cell == Constants.CELL_X)           
                        possibleWinX++;

                   if(cell == Constants.CELL_O)
                        possibleWin0++;

                } // end inner for loop

                if (possibleWinX == BOARD_SIZE)
                    nrOrthWinnableX++;
                if (possibleWin0 == BOARD_SIZE)
                    nrOrthWinnableO++;

        } 
    } 

      // Process Diagonals for layers-plane. i.e  layers are fixed
    // for each iteration. Going from coordinates (0,0,layer) to (3,3,layer) for each layer. 
       for (int layer=0; layer<BOARD_SIZE; layer++) { 

           int possibleWinX =0; int possibleWin0 =0;
           
             for (int row=0, col=0; row < BOARD_SIZE; row++,col++) {

                   int cell = gameState.at(row,col,layer);

                   if (cell == Constants.CELL_EMPTY){
                        possibleWinX++;
                        possibleWin0++;
                   }

                   if (cell == Constants.CELL_X)           
                        possibleWinX++;

                   if(cell == Constants.CELL_O)
                        possibleWin0++;

            } // end outer for loop   
         if (possibleWinX == BOARD_SIZE)
                zDiagX++;
         if (possibleWin0 == BOARD_SIZE)
                zDiagO++;

        }
      
    // Process Diagonals for layers-plane. i.e  layers are fixed
    // for each iteration. Going from (0,3,layer) to (3,0,layer) for each layer. 
       for (int layer=0; layer<BOARD_SIZE; layer++) { 

            int possibleWinX =0; int possibleWin0 =0;

             for (int row=0, col=BOARD_SIZE-1; row < BOARD_SIZE;  row++,col--) {

                   int cell = gameState.at(row,col, layer);

                   if (cell == Constants.CELL_EMPTY){
                        possibleWinX++;
                        possibleWin0++;
                   }

                   if (cell == Constants.CELL_X)           
                        possibleWinX++;

                   if(cell == Constants.CELL_O)
                        possibleWin0++;  

            } // end outer for loop 
            if (possibleWinX == BOARD_SIZE)
                zDiagX++;
            if (possibleWin0 == BOARD_SIZE)
                zDiagO++;
        }

    // Process Diagonals for column-plane. i.e  columns are fixed
    // for each iteration. Going from (0,col,0) to (3,col,3) 
       for (int col=0; col<BOARD_SIZE; col++) { 

           int possibleWinX =0; int possibleWin0 =0;
           
             for (int row=0, layer=0; row < BOARD_SIZE; row++,layer++) {

                   int cell = gameState.at(row,col,layer);

                   if (cell == Constants.CELL_EMPTY){
                        possibleWinX++;
                        possibleWin0++;
                   }

                   if (cell == Constants.CELL_X)           
                        possibleWinX++;

                   if(cell == Constants.CELL_O)
                        possibleWin0++;

            } // end outer for loop   
         if (possibleWinX == BOARD_SIZE)
                yDiagX++;
         if (possibleWin0 == BOARD_SIZE)
                yDiagO++;

        }   

     // Process  Diagonals for column-plane. i.e  columns are fixed
    // for each iteration. Going from (0,col,3) to (3,col,0)
       for (int col=0; col<BOARD_SIZE; col++) { 

            int possibleWinX =0; int possibleWin0 =0;

             for (int row=0, layer=BOARD_SIZE-1; row < BOARD_SIZE;  row++,layer--) {

                   int cell = gameState.at(row,col, layer);

                   if (cell == Constants.CELL_EMPTY){
                        possibleWinX++;
                        possibleWin0++;
                   }

                   if (cell == Constants.CELL_X)           
                        possibleWinX++;

                   if(cell == Constants.CELL_O)
                        possibleWin0++;  

            } // end outer for loop 
            if (possibleWinX == BOARD_SIZE)
                yDiagX++;
            if (possibleWin0 == BOARD_SIZE)
                yDiagO++;
        }

    // Process Diagonals for row-plane. i.e  rows are fixed
    // for each iteration. Going from coordinates (row,0,0) to (row,3,3)
       for (int row=0; row<BOARD_SIZE; row++) { 

           int possibleWinX =0; int possibleWin0 =0;
           
             for (int col=0, layer=0; col < BOARD_SIZE; col++,layer++) {

                   int cell = gameState.at(row,col,layer);

                   if (cell == Constants.CELL_EMPTY){
                        possibleWinX++;
                        possibleWin0++;
                   }

                   if (cell == Constants.CELL_X)           
                        possibleWinX++;

                   if(cell == Constants.CELL_O)
                        possibleWin0++;

            } // end outer for loop   
         if (possibleWinX == BOARD_SIZE)
                xDiagX++;
         if (possibleWin0 == BOARD_SIZE)
                xDiagO++;

        }   

     // Process  Diagonals for row-plane. i.e  rows are fixed
    // for each iteration. Going from coordinates (row,0,3) to (row,3,0)
       for (int row=0; row<BOARD_SIZE; row++) { 

            int possibleWinX =0; int possibleWin0 =0;

             for (int col=0, layer=BOARD_SIZE-1; col < BOARD_SIZE;  col++,layer--) {

                   int cell = gameState.at(row,col, layer);

                   if (cell == Constants.CELL_EMPTY){
                        possibleWinX++;
                        possibleWin0++;
                   }

                   if (cell == Constants.CELL_X)           
                        possibleWinX++;

                   if(cell == Constants.CELL_O)
                        possibleWin0++;  

            } // end outer for loop 
            if (possibleWinX == BOARD_SIZE)
                xDiagX++;
            if (possibleWin0 == BOARD_SIZE)
                xDiagO++;
        }

    // Process the FOUR MAIN Diagonals 
    for (int diag =0; diag < BOARD_SIZE; diag++) {

        int possibleWinX =0; int possibleWin0 =0;
        for (int pos=0; pos < BOARD_SIZE; pos++) {

           int cell = gameState.at(MDIAGS[diag][pos][0], MDIAGS[diag][pos][1],MDIAGS[diag][pos][2]);

           if (cell == Constants.CELL_EMPTY){
                possibleWinX++;
                possibleWin0++;
           }

           if (cell == Constants.CELL_X)           
                possibleWinX++;

           if(cell == Constants.CELL_O)
                possibleWin0++;  

        } // end for loop 
    
        if (possibleWinX == BOARD_SIZE)
            nrDiagWinnableX++;
        if (possibleWin0 == BOARD_SIZE)
            nrDiagWinnableO++;
    }

    int xWins = nrOrthWinnableX + nrRowsWinnableX + nrColWinnableX + nrDiagWinnableX +zDiagX+yDiagX+xDiagX;     ;
    int oWins = nrOrthWinnableO + nrRowsWinnableO + nrColWinnableO + nrDiagWinnableO +zDiagO+yDiagO+xDiagO;
    int res = xWins - oWins;


  /*  System.err.println("Reached evaluate state with row X wins = " + nrRowsWinnableX ); 
    System.err.println("Reached evaluate state with row O wins = " + nrRowsWinnableO );

    System.err.println("Reached evaluate state with column X wins = " + nrColWinnableX ); 
    System.err.println("Reached evaluate state with column O wins = " + nrColWinnableO );

    System.err.println("Reached evaluate state with Main Diag X wins = " + nrDiagWinnableX ); 
    System.err.println("Reached evaluate state with Main Diag O wins = " + nrDiagWinnableO );
*/
  /*  System.err.println("Reached evaluate state with Ortho X wins = " + nrOrthWinnableX ); 
    System.err.println("Reached evaluate state with Ortho O wins = " + nrOrthWinnableO );

    System.err.println("Reached evaluate state with zplane Diag X wins = " + zDiagX ); 
    System.err.println("Reached evaluate state with zplane Diag O wins = " + zDiagO );

    System.err.println("Reached evaluate state with yplane Diag X wins = " + yDiagX ); 
    System.err.println("Reached evaluate state with yplane Diag O wins = " + yDiagO );
    
    System.err.println("Reached evaluate state with xplane Diag X wins = " + xDiagX ); 
    System.err.println("Reached evaluate state with xplane Diag O wins = " + xDiagO );
*/
    return res;

    } 

    private int evalTerminalState(final GameState gState)
    {

        if (gState.isXWin()) {
            return POSINFINITY;
        }
        else if (gState.isOWin()) {
            return NEGINFINITY;
        }
        else {
            return 0;   // terminalState.isDraw()
        }
          
    }
}
