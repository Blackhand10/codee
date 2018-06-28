package Nocitv5;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.IOException;
import java.util.concurrent.ExecutionException;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import javax.swing.JProgressBar;
import javax.swing.filechooser.FileNameExtensionFilter;
import Nocitv5.NocitAlgorithm;
import Nocitv5.NocitInputs;
import Nocitv5.NocitWindow;
import org.apache.commons.math.FunctionEvaluationException;
import org.apache.commons.math.optimization.OptimizationException;

public class NocitControl implements ActionListener {
   private NocitWindow ui;
   private NocitInputs inputs;
   private NocitControl.NocitThread nocitThread;
   private JProgressBar progressBar;
   private JFileChooser fileChooser = new JFileChooser();
   private NocitAlgorithm nocitAlgorithm;
   private final FileNameExtensionFilter textFileFilter = new FileNameExtensionFilter("txt", new String[]{"txt"});
   private final FileNameExtensionFilter csvFileFilter = new FileNameExtensionFilter("csv", new String[]{"csv"});

   public NocitControl(NocitWindow ui, NocitInputs inputs) {
      this.ui = ui;
      this.inputs = inputs;
      this.nocitThread = new NocitControl.NocitThread();
      this.nocitAlgorithm = new NocitAlgorithm();
      this.progressBar = ui.getProgressBar();
      ui.setVisible(true);
      ui.addActionListenerToAllButtons(this);
   }

   public void runNocit() {
      this.progressBar.setIndeterminate(true);
      this.progressBar.setString("NOCIt is running");
      boolean stoppedBeforeCompletion = false;

      try {
         this.nocitAlgorithm.runNocit(this.inputs.getSampleFile(), this.inputs.getFrequencyFile(), this.inputs.getOutputFile(), this.inputs.getCalibrationFile(), this.inputs.getMaxNumberOfContributors(), this.inputs.getDNAMass(), this.progressBar);
         this.nocitAlgorithm.printNocitResults();
      } catch (IOException | FunctionEvaluationException | OptimizationException | ExecutionException var3) {
         System.out.println(var3.getMessage());
      } catch (InterruptedException var4) {
         stoppedBeforeCompletion = true;
         this.nocitAlgorithm.printInterruptedNocitResults();
      }

      this.nocitThread = new NocitControl.NocitThread();
      if(stoppedBeforeCompletion) {
         JOptionPane.showMessageDialog(this.ui, "The program was stopped before completing computation");
      } else {
         JOptionPane.showMessageDialog(this.ui, "The program has finished running");
      }

      this.ui.setStartButtonEnabled(Boolean.TRUE);
      this.progressBar.setString("");
      this.progressBar.setIndeterminate(false);
   }

   public void actionPerformed(ActionEvent e) {
      String var2 = e.getActionCommand();
      byte var3 = -1;
      switch(var2.hashCode()) {
      case -1005512447:
         if(var2.equals("output")) {
            var3 = 3;
         }
         break;
      case -909675094:
         if(var2.equals("sample")) {
            var3 = 2;
         }
         break;
      case -70023844:
         if(var2.equals("frequency")) {
            var3 = 0;
         }
         break;
      case 109757538:
         if(var2.equals("start")) {
            var3 = 4;
         }
         break;
      case 1421318634:
         if(var2.equals("calibration")) {
            var3 = 1;
         }
      }

      File start1;
      switch(var3) {
      case 0:
         start1 = this.getFile("Select Frequency File", this.csvFileFilter);
         if(start1 != null) {
            this.ui.setFrequencyFileSelectionText(start1.getName());
            this.inputs.setFrequencyFile(start1);
         }
         break;
      case 1:
         start1 = this.getFile("Select Calibration File", this.csvFileFilter);
         if(start1 != null) {
            this.ui.setCalibrationFileSelectionText(start1.getName());
            this.inputs.setCalibrationFile(start1);
         }
         break;
      case 2:
         start1 = this.getFile("Select Sample File", this.csvFileFilter);
         if(start1 != null) {
            this.ui.setSampleFileSelectionText(start1.getName());
            this.inputs.setSampleFile(start1);
         }
         break;
      case 3:
         start1 = this.makeFile("Select/Create Output File", this.textFileFilter);
         if(start1 != null) {
            this.ui.setOutputFileSelectionText(start1.getName());
            this.inputs.setOutputFile(start1);
         }
         break;
      case 4:
         boolean start = true;
         if(this.inputs.getOutputFile() != null && this.inputs.getOutputFile().exists() && this.inputs.getOutputFile().length() != 0L) {
            int paneReturn = JOptionPane.showConfirmDialog(this.ui, "The file " + this.inputs.getOutputFile().getName() + " contains NOCIt program output. Do you want to overwrite it?", "Overwrite file?", 0);
            if(paneReturn != 0) {
               start = false;
            }
         }

         if(start && this.inputs.updateNumericFields(this.ui.getDNAMassSelectionText(), this.ui.getMaxNoContibutorsText()) && this.inputs.allFieldsValid()) {
            if(this.inputs.getOutputFile() == null) {
               this.inputs.setOutputFile(new File(this.getAbsoluteFilePathWithoutExtension(this.inputs.getSampleFile()) + "_NOCIt_output.txt"));
               this.ui.setOutputFileSelectionText(this.inputs.getOutputFile().getName());
            }

            this.nocitThread.start();
            this.ui.setStartButtonEnabled(Boolean.FALSE);
         }

         this.displayNumericFieldErrors();
         break;
      default:
         this.nocitThread.interrupt();
         this.progressBar.setString("");
         this.progressBar.setIndeterminate(false);
      }

   }

   private File getFile(String message, FileNameExtensionFilter filter) {
      this.fileChooser.setDialogTitle(message);
      this.fileChooser.setFileFilter(filter);
      int chooserReturn = this.fileChooser.showOpenDialog(this.ui);
      File selectedFile = null;
      if(chooserReturn == 0) {
         selectedFile = this.fileChooser.getSelectedFile();
      }

      return selectedFile;
   }

   private File makeFile(String message, FileNameExtensionFilter filter) {
      this.fileChooser.setFileFilter(filter);
      this.fileChooser.setDialogTitle(message);
      int chooserReturn = this.fileChooser.showSaveDialog(this.ui);
      File createdFile = null;
      if(chooserReturn == 0) {
         createdFile = this.fileChooser.getSelectedFile();
         if(createdFile.exists()) {
            int pathName = JOptionPane.showConfirmDialog(this.ui, "The file \"" + createdFile.getName() + "\" exists and will be overwritten by the output of the program. Do you want to continue?", "Overwrite file?", 0);
            if(pathName == 0) {
               createdFile.delete();

               try {
                  createdFile.createNewFile();
               } catch (IOException var9) {
                  System.out.println(var9.getMessage());
               }
            } else {
               createdFile = null;
            }
         } else {
            String pathName1 = createdFile.getAbsolutePath();
            int begin = pathName1.length() - 4;
            if(!createdFile.getAbsolutePath().substring(begin).equals(".txt")) {
               createdFile = new File(createdFile.getAbsolutePath().concat(".txt"));
            }

            try {
               createdFile.createNewFile();
            } catch (IOException var8) {
               System.out.println(var8.getMessage());
            }
         }
      }

      return createdFile;
   }

   private String getAbsoluteFilePathWithoutExtension(File file) {
      int end = file.getAbsolutePath().length() - 4;
      return this.inputs.getSampleFile().getAbsolutePath().substring(0, end);
   }

   private void displayNumericFieldErrors() {
      String dnaMassError = this.inputs.getDnaMassError();
      String maxNumberOfContributorsError = this.inputs.getMaxNumberOfContributorsError();
      if(!dnaMassError.equals("") && !maxNumberOfContributorsError.equals("")) {
         JOptionPane.showMessageDialog(this.ui, "Inputs for \"DNA mass\" and \"maximum number of contributors\" are not valid");
      } else if(!dnaMassError.equals("") && maxNumberOfContributorsError.equals("")) {
         JOptionPane.showMessageDialog(this.ui, dnaMassError);
      } else if(dnaMassError.equals("") && !maxNumberOfContributorsError.equals("")) {
         JOptionPane.showMessageDialog(this.ui, maxNumberOfContributorsError);
      }

      this.inputs.setDnaMassError("");
      this.inputs.setMaxNumberOfContributorsError("");
   }

   class NocitThread extends Thread {
      public void run() {
         NocitControl.this.runNocit();
      }
   }
}
