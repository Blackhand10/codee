package Nocitv5;

import java.awt.Component;
import java.io.File;
import javax.swing.JOptionPane;

public class NocitInputs {
   private File calibrationFile = null;
   private File frequencyFile = null;
   private File sampleFile = null;
   private File outputFile = null;
   private int maxNumberOfContributors = 0;
   private double dnaMass = 0.0D;
   private String dnaMassError = "";
   private String maxNumberOfContributorsError = "";

   public boolean updateNumericFields(String mass, String maxNumberContributors) {
      boolean fieldsValid = true;

      try {
         this.dnaMass = Double.parseDouble(mass);
      } catch (IllegalArgumentException var6) {
         this.dnaMassError = "Invalid DNA mass value.";
         fieldsValid = false;
      }

      try {
         this.maxNumberOfContributors = Integer.parseInt(maxNumberContributors);
      } catch (IllegalArgumentException var5) {
         this.maxNumberOfContributorsError = "Invalid value for maximum number of contributors. ";
         fieldsValid = false;
      }

      return fieldsValid;
   }

   public boolean allFieldsValid() {
      if(this.calibrationFile == null) {
         JOptionPane.showMessageDialog((Component)null, "A calibration file was not provided");
         return false;
      } else if(this.frequencyFile == null) {
         JOptionPane.showMessageDialog((Component)null, "A frequency file was not provided");
         return false;
      } else if(this.sampleFile == null) {
         JOptionPane.showMessageDialog((Component)null, "A sample file was not provided");
         return false;
      } else if(this.dnaMass < 0.0D) {
         this.dnaMassError = "Invalid DNA mass value.";
         return false;
      } else if(this.maxNumberOfContributors >= 0 && this.maxNumberOfContributors <= 5) {
         return true;
      } else {
         this.maxNumberOfContributorsError = "Invalid value for maximum number of contributors.";
         return false;
      }
   }

   public void setDnaMassError(String error) {
      this.dnaMassError = error;
   }

   public void setMaxNumberOfContributorsError(String error) {
      this.maxNumberOfContributorsError = error;
   }

   public void setCalibrationFile(File calibrationFile) {
      this.calibrationFile = calibrationFile;
   }

   public void setFrequencyFile(File frequencyFile) {
      this.frequencyFile = frequencyFile;
   }

   public void setSampleFile(File sampleFile) {
      this.sampleFile = sampleFile;
   }

   public void setOutputFile(File outputFile) {
      this.outputFile = outputFile;
   }

   public void setMaxNumberOfContributors(int maxNumberOfContributors) {
      this.maxNumberOfContributors = maxNumberOfContributors;
   }

   public void setDNAMass(double mass) {
      this.dnaMass = mass;
   }

   public File getCalibrationFile() {
      return this.calibrationFile;
   }

   public File getFrequencyFile() {
      return this.frequencyFile;
   }

   public File getSampleFile() {
      return this.sampleFile;
   }

   public File getOutputFile() {
      return this.outputFile;
   }

   public int getMaxNumberOfContributors() {
      return this.maxNumberOfContributors;
   }

   public double getDNAMass() {
      return this.dnaMass;
   }

   public String getDnaMassError() {
      return this.dnaMassError;
   }

   public String getMaxNumberOfContributorsError() {
      return this.maxNumberOfContributorsError;
   }
}
