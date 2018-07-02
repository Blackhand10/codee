/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
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
import org.apache.commons.math.FunctionEvaluationException;
import org.apache.commons.math.optimization.OptimizationException;
import java.awt.Component;
import java.io.File;
import javax.swing.JOptionPane;


/**
 *
 * @author emada
 */
public class running {
    public static void main(String[] args){
    NocitAlgorithm no = new NocitAlgorithm();
    File calibrationFile = new File("C:/Users/emada/Documents/NOCIT/NOCItpackage - v15/NOCIt Tutorial 1 v1.5.1/NOCIt Tutorial 1 - Calibration file_10sec.csv");
    File frequencyFile = new File("C:/Users/emada/Documents/NOCIT/NOCItpackage - v15/NOCIt Tutorial 1 v1.5.1/NOCIt Tutorial 1 - Caucasian IDPlus Freq file.csv");
    File sampleFile = new File("C:/Users/emada/Documents/NOCIT/NOCItpackage - v15/NOCIt Tutorial 1 v1.5.1/NOCIt Tutorial 1 - SingleSource-0.25IP-10sec.csv");
    File outputFile = new File("C:/Users/emada/Documents/NOCIT/NOCItpackage - v15/NOCIt Tutorial 1 v1.5.1/NOCIt Tutorial 1 v1.5.1 - Answers/output.txt");
    int maxNumberOfContributors = 1;
    double dnaMass = .25;
    try {
         no.runNocit(sampleFile, frequencyFile, outputFile, calibrationFile, maxNumberOfContributors, dnaMass);
         no.printNocitResults();

      } catch (IOException | FunctionEvaluationException | OptimizationException | ExecutionException var3) {
         System.out.println(var3.getMessage());
      } catch (InterruptedException var4) {
         no.printInterruptedNocitResults();
      }
        }
    
}
