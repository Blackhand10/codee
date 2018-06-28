package Nocitv5;

import Nocitv5.NocitControl;
import Nocitv5.NocitInputs;
import Nocitv5.NocitWindow;

public class NocitMain {
   public static void main(String[] args) {
      NocitWindow ui = new NocitWindow();
      NocitInputs inputs = new NocitInputs();
      new NocitControl(ui, inputs);
   }
}
