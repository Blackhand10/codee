package Nocitv5;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.Iterator;
import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JTextField;
import javax.swing.border.EmptyBorder;
import javax.swing.border.EtchedBorder;

public class NocitWindow extends JFrame {
   private JPanel basePanel;
   private JPanel filesPanel;
   private JPanel numericInputPanel;
   private JPanel DNAMassPanel;
   private JPanel maxNumberContibutorsPanel;
   private JPanel buttonsPanel;
   private JPanel panelForFilesPanelBorder;
   private JPanel progressBarPanel;
   private JLabel calibrationFileLabel;
   private JLabel sampleFileLabel;
   private JLabel frequencyFileLabel;
   private JLabel outputFileLabel;
   private JLabel DNAMassLabel;
   private JLabel maxNumberContributorsLabel;
   private JTextField DNAMassSelection;
   private JTextField maxNumberContributorsSelection;
   private final Dimension fileSelectionValueLabelDimension = new Dimension(200, 20);
   private JLabel sampleFileSelectionValueLabel;
   private JLabel frequencyFileSelectionValueLabel;
   private JLabel calibrationFileSelectionValueLabel;
   private JLabel outputFileSelectionValueLabel;
   private ArrayList buttons;
   private JProgressBar progressBar;
   private JButton calibrationFileSelectionButton;
   private JButton sampleFileSelectionButton;
   private JButton frequencyFileSelectionButton;
   private JButton createOutputFileButton;
   private JButton startButton;
   private JButton stopButton;

   public NocitWindow() {
      this.buttons = new ArrayList();
      this.setupPanels();
      this.setupFrame();
   }

   public NocitWindow(ActionListener listener) {
      this.setupPanels();
      this.setupFrame();
      this.calibrationFileSelectionButton.addActionListener(listener);
   }

   public JButton getSampleFileSelectionButton() {
      return this.sampleFileSelectionButton;
   }

   private void setupPanels() {
      this.basePanel = new JPanel();
      this.basePanel.setOpaque(true);
      this.basePanel.setLayout(new BoxLayout(this.basePanel, 3));
      this.basePanel.setBorder(new EmptyBorder(5, 5, 5, 5));
      this.calibrationFileLabel = new JLabel("Calibration file: ");
      this.sampleFileLabel = new JLabel("Sample file:  ");
      this.outputFileLabel = new JLabel("Output file:");
      this.frequencyFileLabel = new JLabel("Frequency file:");
      this.DNAMassLabel = new JLabel("Sample DNA input:");
      this.maxNumberContributorsLabel = new JLabel("Maximum number of contributors (0-5):");
      this.calibrationFileSelectionValueLabel = new JLabel(" ", 0);
      this.calibrationFileSelectionValueLabel.setName("calibrationFileSelection");
      this.calibrationFileSelectionValueLabel.setPreferredSize(this.fileSelectionValueLabelDimension);
      this.calibrationFileSelectionValueLabel.setBorder(new EtchedBorder());
      this.calibrationFileSelectionValueLabel.setOpaque(true);
      this.calibrationFileSelectionValueLabel.setBackground(Color.WHITE);
      this.sampleFileSelectionValueLabel = new JLabel("", 0);
      this.sampleFileSelectionValueLabel.setName("sampleFileSelection");
      this.sampleFileSelectionValueLabel.setPreferredSize(this.fileSelectionValueLabelDimension);
      this.sampleFileSelectionValueLabel.setBorder(new EtchedBorder());
      this.sampleFileSelectionValueLabel.setOpaque(true);
      this.sampleFileSelectionValueLabel.setBackground(Color.WHITE);
      this.frequencyFileSelectionValueLabel = new JLabel("", 0);
      this.frequencyFileSelectionValueLabel.setName("frequencyFileSelection");
      this.frequencyFileSelectionValueLabel.setPreferredSize(this.fileSelectionValueLabelDimension);
      this.frequencyFileSelectionValueLabel.setBorder(new EtchedBorder());
      this.frequencyFileSelectionValueLabel.setOpaque(true);
      this.frequencyFileSelectionValueLabel.setBackground(Color.WHITE);
      this.outputFileSelectionValueLabel = new JLabel("", 0);
      this.outputFileSelectionValueLabel.setName("outputFileSelection");
      this.outputFileSelectionValueLabel.setPreferredSize(this.fileSelectionValueLabelDimension);
      this.outputFileSelectionValueLabel.setBorder(new EtchedBorder());
      this.outputFileSelectionValueLabel.setOpaque(true);
      this.outputFileSelectionValueLabel.setBackground(Color.WHITE);
      this.DNAMassSelection = new JTextField("", 5);
      this.DNAMassSelection.setName("DNAMassSelection");
      this.maxNumberContributorsSelection = new JTextField("", 5);
      this.maxNumberContributorsSelection.setName("maxNumberContributorsSelection");
      this.calibrationFileSelectionButton = new JButton("Browse");
      this.calibrationFileSelectionButton.setActionCommand("calibration");
      this.buttons.add(this.calibrationFileSelectionButton);
      this.sampleFileSelectionButton = new JButton("Browse");
      this.sampleFileSelectionButton.setActionCommand("sample");
      this.buttons.add(this.sampleFileSelectionButton);
      this.frequencyFileSelectionButton = new JButton("Browse");
      this.frequencyFileSelectionButton.setActionCommand("frequency");
      this.buttons.add(this.frequencyFileSelectionButton);
      this.createOutputFileButton = new JButton("Browse");
      this.createOutputFileButton.setActionCommand("output");
      this.buttons.add(this.createOutputFileButton);
      this.startButton = new JButton("START");
      this.startButton.setActionCommand("start");
      this.buttons.add(this.startButton);
      this.stopButton = new JButton("STOP");
      this.stopButton.setActionCommand("stop");
      this.buttons.add(this.stopButton);
      this.DNAMassPanel = new JPanel();
      this.DNAMassPanel.setLayout(new FlowLayout(0));
      this.DNAMassPanel.add(this.DNAMassLabel);
      this.DNAMassPanel.add(this.DNAMassSelection);
      this.maxNumberContibutorsPanel = new JPanel();
      this.maxNumberContibutorsPanel.setLayout(new FlowLayout(0));
      this.maxNumberContibutorsPanel.add(this.maxNumberContributorsLabel);
      this.maxNumberContibutorsPanel.add(this.maxNumberContributorsSelection);
      this.numericInputPanel = new JPanel();
      this.numericInputPanel.setLayout(new BoxLayout(this.numericInputPanel, 3));
      this.numericInputPanel.setBorder(BorderFactory.createCompoundBorder(new EmptyBorder(10, 10, 5, 10), BorderFactory.createEtchedBorder()));
      this.numericInputPanel.add(this.DNAMassPanel);
      this.numericInputPanel.add(this.maxNumberContibutorsPanel);
      this.filesPanel = new JPanel(new GridBagLayout());
      this.filesPanel.setBorder(BorderFactory.createEmptyBorder(3, 3, 3, 3));
      GridBagConstraints c = new GridBagConstraints();
      c.gridx = 0;
      c.gridy = 0;
      c.insets = new Insets(3, 3, 3, 3);
      c.fill = 0;
      this.filesPanel.add(this.calibrationFileLabel, c);
      c.gridy = 1;
      this.filesPanel.add(this.frequencyFileLabel, c);
      c.gridy = 2;
      this.filesPanel.add(this.sampleFileLabel, c);
      c.gridy = 3;
      this.filesPanel.add(this.outputFileLabel, c);
      c.gridx = 1;
      c.gridy = 0;
      c.fill = 1;
      this.filesPanel.add(this.calibrationFileSelectionValueLabel, c);
      c.gridy = 1;
      this.filesPanel.add(this.frequencyFileSelectionValueLabel, c);
      c.gridy = 2;
      this.filesPanel.add(this.sampleFileSelectionValueLabel, c);
      c.gridy = 3;
      this.filesPanel.add(this.outputFileSelectionValueLabel, c);
      c.gridx = 2;
      c.gridy = 0;
      c.fill = 0;
      this.filesPanel.add(this.calibrationFileSelectionButton, c);
      c.gridy = 1;
      this.filesPanel.add(this.frequencyFileSelectionButton, c);
      c.gridy = 2;
      this.filesPanel.add(this.sampleFileSelectionButton, c);
      c.gridy = 3;
      this.filesPanel.add(this.createOutputFileButton, c);
      this.panelForFilesPanelBorder = new JPanel();
      this.panelForFilesPanelBorder.setBorder(BorderFactory.createCompoundBorder(new EmptyBorder(5, 10, 5, 10), BorderFactory.createEtchedBorder()));
      this.panelForFilesPanelBorder.add(this.filesPanel);
      this.progressBar = new JProgressBar();
      this.progressBar.setStringPainted(true);
      this.progressBar.setString("");
      this.progressBarPanel = new JPanel();
      this.progressBarPanel.add(this.progressBar);
      this.progressBarPanel.setBorder(BorderFactory.createCompoundBorder(new EmptyBorder(5, 10, 5, 10), BorderFactory.createEtchedBorder()));
      this.buttonsPanel = new JPanel();
      this.buttonsPanel.setLayout(new FlowLayout(1));
      this.buttonsPanel.setBorder(BorderFactory.createCompoundBorder(new EmptyBorder(5, 10, 10, 10), BorderFactory.createEtchedBorder()));
      this.buttonsPanel.add(this.startButton);
      this.buttonsPanel.add(this.stopButton);
      this.basePanel.add(this.panelForFilesPanelBorder);
      this.basePanel.add(this.numericInputPanel);
      this.basePanel.add(this.progressBarPanel);
      this.basePanel.add(this.buttonsPanel);
   }

   private void setupFrame() {
      this.setContentPane(this.basePanel);
      this.setDefaultCloseOperation(3);
      this.setResizable(false);
      this.setTitle("NOCIt");
      this.pack();
   }

   public void addActionListenerToAllButtons(ActionListener listener) {
      Iterator i$ = this.buttons.iterator();

      while(i$.hasNext()) {
         JButton b = (JButton)i$.next();
         b.addActionListener(listener);
      }

   }

   public void setCalibrationFileSelectionText(String selection) {
      this.calibrationFileSelectionValueLabel.setText(selection);
   }

   public String getCalibrationFileSelectionText() {
      return this.calibrationFileSelectionValueLabel.getText();
   }

   public void setFrequencyFileSelectionText(String selection) {
      this.frequencyFileSelectionValueLabel.setText(selection);
   }

   public String getFrequencyFileSelectionText() {
      return this.frequencyFileSelectionValueLabel.getText();
   }

   public void setSampleFileSelectionText(String selection) {
      this.sampleFileSelectionValueLabel.setText(selection);
   }

   public String getSampleFileSelectionText() {
      return this.sampleFileSelectionValueLabel.getText();
   }

   public JProgressBar getProgressBar() {
      return this.progressBar;
   }

   public void setOutputFileSelectionText(String selection) {
      this.outputFileSelectionValueLabel.setText(selection);
   }

   public String getOutputFileSelectionText() {
      return this.outputFileSelectionValueLabel.getText();
   }

   public String getDNAMassSelectionText() {
      return this.DNAMassSelection.getText();
   }

   public String getMaxNoContibutorsText() {
      return this.maxNumberContributorsSelection.getText();
   }

   public void setStartButtonEnabled(Boolean b) {
      this.startButton.setEnabled(b.booleanValue());
   }
}
