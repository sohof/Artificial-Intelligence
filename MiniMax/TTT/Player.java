import java.util.*;

public class Player {

    private static final int NEGINFINITY = Integer.MIN_VALUE;
    private static final int POSINFINITY = Integer.MAX_VALUE;
    private long TIMELIMIT = 850000000; //timelmit for each board is max of 0.25 seconds
    private int depth=1;
    /**
     * Performs a move
     * @param gameState the current state of the board
     * @param deadline time before which we must have returned
     * @return the next state the board is in after our move
     */

    public GameState play(final GameState gameState, final Deadline deadline) {
        Vector<GameState> nextStates = new Vector<GameState>();

        gameState.findPossibleMoves(nextStates);
        int vSize = nextStates.size();

        if (nextStates.size() == 0) {
            // Must play "pass" move if there are no other moves possible.
            return new GameState(gameState, new Move());
        }

        int argmax=-1;      
        int bestValue = NEGINFINITY;
        int maxDepth =3;
        for (int i =0; i < vSize; i++){
               // System.err.println(nextStates.get(i).toString(Constants.CELL_X));
                int alpha =  Integer.MIN_VALUE;
                int beta =  Integer.MAX_VALUE;
                int res = minValue(nextStates.get(i),1,maxDepth, alpha, beta, deadline);
                //System.err.println("minValue for given state " + res);

                if (res > bestValue){
                    bestValue = res;
                    argmax = i;
                }

            
            //System.err.println("Time remaining " + deadline.timeUntil() + " Reached depth " + maxDepth);
        }
        return nextStates.get(argmax);
     
    }
    private int minValue(GameState gState, int depth, int maxDepth, int alpha, int beta, Deadline deadline) {

         if (gState.isEOG()) {

            return evalTerminalState(gState);
         }

         if (depth == maxDepth) {
            return evalState(gState, depth);
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
            return evalState(gState,depth);
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
    private int heur(int possibleWin){

        if (possibleWin == 1)
            return 1;
        if (possibleWin == 2)
            return 10;
        if (possibleWin == 3)
            return 100; 

        System.err.println("heur function error"); 
        return 0;
    }


    private int evalState(final GameState gameState, int depth) 
    {

        System.err.println("Eval state reached depth " + depth); 

        final int BOARD_SIZE = gameState.BOARD_SIZE;
        int nrRowsWinnableX =0; 
        int nrRowsWinnableO =0;   
        int nrColWinnableX =0; 
        int nrColWinnableO =0; 
        int nrDiagWinnableX =0; 
        int nrDiagWinnableO =0;

        // Process all the rows
        for (int row=0; row < BOARD_SIZE; row++) {

            int pWinX =0; int pWinO =0;

            for (int col=0; col < BOARD_SIZE; col++) {

               int cell = gameState.at(row,col);

               if (cell == Constants.CELL_X)           
                    pWinX++;

               if(cell == Constants.CELL_O)
                    pWinO++;

               if (pWinO >0 && pWinX>0)
                     break;   
            } // end inner for loop
                if (pWinX > 0)
                    nrRowsWinnableX += heur(pWinX);
                
                if (pWinO > 0)
                    nrRowsWinnableO += heur(pWinO);

        } // end outer for loop

        // Process the columns
        for (int col =0; col< BOARD_SIZE; col++) {

            int pWinX =0; int pWinO =0;

            for (int row=0; row < BOARD_SIZE; row++) {

               int cell = gameState.at(row,col);

               if (cell == Constants.CELL_X)           
                    pWinX++;

               if(cell == Constants.CELL_O)
                    pWinO++;

            } // end inner for loop

            if (pWinX > 0)
                    nrColWinnableX += heur(pWinX);
                
            if (pWinO > 0)
                nrColWinnableO += heur(pWinO);
        } // end outer for loop


      // Process the Left Diagonal
        boolean run = true;
        while (run) {

        int pWinX =0; int pWinO =0;
           
             for (int row=0, col=0; row < BOARD_SIZE; row++,col++) {

                   int cell = gameState.at(row,col);

                   if (cell == Constants.CELL_X)           
                        pWinX++;

                   if(cell == Constants.CELL_O)
                        pWinO++;

            } // end outer for loop   
            if (pWinX > 0)
                    nrDiagWinnableX += heur(pWinX);
                
            if (pWinO > 0)
                    nrDiagWinnableO += heur(pWinO);

        run = false;
        }

      // Process the Right Diagonal
        boolean run2 =true;
        while (run2) {

            int pWinX =0; int pWinO =0;

             for (int row=0, col=BOARD_SIZE-1; row < BOARD_SIZE;  row++,col--) {

                   int cell = gameState.at(row,col);

                   if (cell == Constants.CELL_X)           
                        pWinX++;

                   if(cell == Constants.CELL_O)
                        pWinO++;  

            } // end outer for loop 
          if (pWinX > 0)
                    nrDiagWinnableX += heur(pWinX);
                
            if (pWinO > 0)
                    nrDiagWinnableO += heur(pWinO);
            run2 = false;
        }

    int res = ((nrRowsWinnableX + nrColWinnableX + nrDiagWinnableX) - (nrRowsWinnableO + nrColWinnableO + nrDiagWinnableO)); 
    //System.err.println("Reached evaluate state with result = " + res ); 
   
    return res;

    } 

    private int evalTerminalState(final GameState gState)
    {

        System.err.println("Reached terminal state"); 

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