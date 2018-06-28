package Nocitv5;
//new comment just to check github sync

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.MathContext;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadLocalRandom;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import javax.swing.JProgressBar;
import org.apache.commons.math.FunctionEvaluationException;
import org.apache.commons.math.optimization.OptimizationException;
import org.apache.commons.math.optimization.fitting.CurveFitter;
import org.apache.commons.math.optimization.fitting.ParametricRealFunction;
import org.apache.commons.math.optimization.general.LevenbergMarquardtOptimizer;
import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;

public class NocitAlgorithm {
   private static HashMap Final_Values = new HashMap();
   private static LinkedHashMap Intervals_Map = new LinkedHashMap();
   private static ConcurrentHashMap Loci_Peaks = new ConcurrentHashMap();
   private static HashMap Loci_Alleles = new HashMap();
   private static HashMap AMEL_Alleles = new HashMap();
   private static ConcurrentHashMap AMEL_Peaks = new ConcurrentHashMap();
   private static HashMap True_Mean_Slope = new HashMap();
   private static HashMap True_Stddev_Slope = new HashMap();
   private static HashMap Noise_Mean_Slope = new HashMap();
   private static HashMap Noise_Stddev_Slope = new HashMap();
   private static HashMap R_Stutter_Mean = new HashMap();
   private static HashMap R_Stutter_Stddev = new HashMap();
   private static HashMap R_Stut_DO = new HashMap();
   private static HashMap F_Stutter_Mean = new HashMap();
   private static HashMap F_Stutter_Stddev = new HashMap();
   private static HashMap F_Stut_DO = new HashMap();
   private static HashMap Locus_DO = new HashMap();
   private static double DNA_Mass;
   private static String[] Sexes;
   private static HashMap Sexes_Genotype;
   private static ArrayList Sample_Loci;
   private String Freq_File;
   private String calibrationPath;
   private String Sample_File_Name;
   private String Output_File_Name;
   private int noc;
   private HashMap n_ans;
   private HashMap n_probdist;
   private HashMap n_locus_value;
   private HashMap n_loci_max_alleles;
   private HashMap n_amel_max_alleles;
   private HashMap time_taken;
   private HashSet Working_Loci;
   private HashSet Loci_RStutter;
   private HashSet Loci_FStutter;
   private HashSet Loci_no_true;
   private HashSet Loci_no_noise;
   private HashSet Loci_not_in_calib_samples;
   private HashSet Loci_not_in_freq_table;

   private static double[] calcSlopeValue(String locusname, double massvalue, HashMap meanslope, HashMap stddevslope) {
      double a_mean = ((double[])meanslope.get(locusname))[0];
      double b_mean = ((double[])meanslope.get(locusname))[1];
      double a_stddev = ((double[])stddevslope.get(locusname))[0];
      double b_stddev = ((double[])stddevslope.get(locusname))[1];
      double mean = a_mean * massvalue + b_mean;
      double stddev = a_stddev * massvalue + b_stddev;
      return new double[]{mean, stddev};
   }

   private static double calcExpValue(String locusname, HashMap locushashmap, double xvalue) {
      double a = ((double[])locushashmap.get(locusname))[0];
      double b = ((double[])locushashmap.get(locusname))[1];
      double expvalue = a * Math.exp(b * xvalue);
      return expvalue;
   }

   private static double calcTwoExpValue(String locusname, HashMap locushashmap, double xvalue) {
      double a = ((double[])locushashmap.get(locusname))[0];
      double b = ((double[])locushashmap.get(locusname))[1];
      double c = ((double[])locushashmap.get(locusname))[2];
      double twoexpvalue = a * Math.exp(b * xvalue) + c;
      return twoexpvalue;
   }

   private static LinkedHashMap intervals(LinkedHashMap freqdict, ArrayList locilist) {
      LinkedHashMap freq_dict = freqdict;
      ArrayList loci_list = locilist;
      LinkedHashMap intervals_map = new LinkedHashMap();
      LinkedHashMap endpoints_map = new LinkedHashMap();

      for(int i$ = 0; i$ < loci_list.size(); ++i$) {
         intervals_map.put(loci_list.get(i$), new LinkedHashMap());
         endpoints_map.put(loci_list.get(i$), new LinkedHashMap());
      }

      Iterator var15 = freqdict.keySet().iterator();

      while(var15.hasNext()) {
         String locus = (String)var15.next();
         int[] keys_array = new int[((LinkedHashMap)freq_dict.get(locus)).size()];
         double[] values_array = new double[((LinkedHashMap)freq_dict.get(locus)).size()];
         int counter = 0;

         int i$1;
         for(Iterator count = ((LinkedHashMap)freq_dict.get(locus)).keySet().iterator(); count.hasNext(); ++counter) {
            i$1 = ((Integer)count.next()).intValue();
            keys_array[counter] = i$1;
            values_array[counter] = ((Double)((LinkedHashMap)freq_dict.get(locus)).get(Integer.valueOf(i$1))).doubleValue();
         }

         for(int var16 = 0; var16 < keys_array.length; ++var16) {
            if(var16 == 0) {
               ((LinkedHashMap)endpoints_map.get(locus)).put(Integer.valueOf(keys_array[var16]), new int[]{0, (int)Math.round(10000.0D * values_array[var16]) - 1});
            } else {
               i$1 = ((int[])((LinkedHashMap)endpoints_map.get(locus)).get(Integer.valueOf(keys_array[var16 - 1])))[1] + 1;
               ((LinkedHashMap)endpoints_map.get(locus)).put(Integer.valueOf(keys_array[var16]), new int[]{i$1, i$1 + (int)Math.round(10000.0D * values_array[var16]) - 1});
            }
         }

         Iterator var17 = ((LinkedHashMap)endpoints_map.get(locus)).keySet().iterator();

         while(var17.hasNext()) {
            int allele = ((Integer)var17.next()).intValue();

            for(int i = ((int[])((LinkedHashMap)endpoints_map.get(locus)).get(Integer.valueOf(allele)))[0]; i <= ((int[])((LinkedHashMap)endpoints_map.get(locus)).get(Integer.valueOf(allele)))[1]; ++i) {
               ((LinkedHashMap)intervals_map.get(locus)).put(Integer.valueOf(i), Integer.valueOf(allele));
            }
         }
      }

      return intervals_map;
   }

   private static double calcLocusPeakHeightsProb(ConcurrentHashMap allpeaks, HashMap mus, HashMap sigmasquareds) {
      double first_term = 1.0D;
      double sum = 0.0D;

      Peak peakobj;
      for(Iterator i$ = allpeaks.values().iterator(); i$.hasNext(); sum += ((double)peakobj.getHeight() - ((Double)mus.get(Integer.valueOf(peakobj.getAllele()))).doubleValue()) * ((double)peakobj.getHeight() - ((Double)mus.get(Integer.valueOf(peakobj.getAllele()))).doubleValue()) / ((Double)sigmasquareds.get(Integer.valueOf(peakobj.getAllele()))).doubleValue()) {
         peakobj = (Peak)i$.next();
         first_term *= 1.0D / Math.sqrt(6.283185307179586D * ((Double)sigmasquareds.get(Integer.valueOf(peakobj.getAllele()))).doubleValue());
      }

      return first_term * Math.exp(-0.5D * sum);
   }

   private static double calcAMELPeakHeightsProb(ConcurrentHashMap allpeaks, HashMap mus, HashMap sigmasquareds) {
      double first_term = 1.0D;
      double sum = 0.0D;

      AMEL_Peak peakobj;
      for(Iterator i$ = allpeaks.values().iterator(); i$.hasNext(); sum += ((double)peakobj.getHeight() - ((Double)mus.get(peakobj.getAllele())).doubleValue()) * ((double)peakobj.getHeight() - ((Double)mus.get(peakobj.getAllele())).doubleValue()) / ((Double)sigmasquareds.get(peakobj.getAllele())).doubleValue()) {
         peakobj = (AMEL_Peak)i$.next();
         first_term *= 1.0D / Math.sqrt(6.283185307179586D * ((Double)sigmasquareds.get(peakobj.getAllele())).doubleValue());
      }

      return first_term * Math.exp(-0.5D * sum);
   }

   private static double[] generateRandomTotal(int numb) {
      int number = numb;
      double[] sorted_array = new double[numb];
      Double[] num_array = new Double[numb - 1];

      for(int count = 1; count <= number - 1; ++count) {
         double i = ThreadLocalRandom.current().nextDouble(1.0D);
         num_array[count - 1] = Double.valueOf(i);
      }

      Arrays.sort(num_array);
      ArrayList num_array1 = new ArrayList(Arrays.asList(num_array));
      ArrayList num_array2 = new ArrayList(Arrays.asList(num_array));
      num_array2.add(Double.valueOf(1.0D));
      num_array1.add(0, Double.valueOf(0.0D));

      for(int var10 = 0; var10 <= number - 1; ++var10) {
         double diff = ((Double)num_array2.get(var10)).doubleValue() - ((Double)num_array1.get(var10)).doubleValue();
         sorted_array[var10] = diff;
      }

      return sorted_array;
   }

   public void printInterruptedNocitResults() {
      BigDecimal check1 = new BigDecimal(0);
      BigDecimal total = new BigDecimal(0);
      int final_n = 0;

      int e;
      for(e = 0; e <= this.noc; ++e) {
         if(this.n_ans.get(Integer.valueOf(e)) == null) {
            final_n = e - 1;
            break;
         }

         total = total.add((BigDecimal)this.n_ans.get(Integer.valueOf(e)));
      }

      if(total.compareTo(check1) == 0) {
         for(e = 0; e <= final_n; ++e) {
            this.n_probdist.put(Integer.valueOf(e), new BigDecimal(0));
         }
      } else {
         for(e = 0; e <= final_n; ++e) {
            BigDecimal output_file = ((BigDecimal)this.n_ans.get(Integer.valueOf(e))).divide(total, MathContext.DECIMAL128);
            if(output_file.compareTo(check1) == 0) {
               this.n_probdist.put(Integer.valueOf(e), new BigDecimal(0));
            } else {
               this.n_probdist.put(Integer.valueOf(e), output_file);
            }
         }
      }

      try {
         FileWriter var27 = new FileWriter(this.Output_File_Name);
         BufferedWriter var26 = new BufferedWriter(var27);
         Throwable var6 = null;

         try {
            var26.write("Calibration File: " + this.calibrationPath + "\n");
            var26.write("Frequency File: " + this.Freq_File + "\n");
            var26.write("Sample File: " + this.Sample_File_Name + "\n");
            var26.write("Sample DNA input: " + DNA_Mass + "ng\n");
            var26.write("\n");
            Iterator x2 = this.Loci_not_in_calib_samples.iterator();

            String i$;
            while(x2.hasNext()) {
               i$ = (String)x2.next();
               var26.write("Locus " + i$ + " was not included in the calculation because it is not present in the calibration samples\n");
            }

            x2 = this.Loci_not_in_freq_table.iterator();

            while(x2.hasNext()) {
               i$ = (String)x2.next();
               var26.write("Locus " + i$ + " was not included in the calculation because it is not present in the frequency table\n");
            }

            x2 = this.Loci_no_true.iterator();

            while(x2.hasNext()) {
               i$ = (String)x2.next();
               var26.write("Locus " + i$ + " was not included in the calculation because it had too few data points for allelic peaks\n");
            }

            x2 = this.Loci_no_noise.iterator();

            while(x2.hasNext()) {
               i$ = (String)x2.next();
               var26.write("Locus " + i$ + " was not included in the calculation because it had too few data points for noise peaks\n");
            }

            x2 = Sample_Loci.iterator();

            while(x2.hasNext()) {
               i$ = (String)x2.next();
               if(!this.Loci_RStutter.contains(i$) && !i$.contentEquals("AMEL")) {
                  var26.write("Locus " + i$ + " had too few data points for reverse stutter \n");
               }
            }

            x2 = Sample_Loci.iterator();

            while(x2.hasNext()) {
               i$ = (String)x2.next();
               if(!this.Loci_FStutter.contains(i$) && !i$.contentEquals("AMEL")) {
                  var26.write("Locus " + i$ + " had too few data points for forward stutter \n");
               }
            }

            var26.write("\n");

            for(int var28 = 0; var28 <= final_n; ++var28) {
               if(this.n_ans.get(Integer.valueOf(var28)) != null) {
                  var26.write("Number of contributors: " + var28 + "\n");
                  var26.write("Time taken: " + this.time_taken.get(Integer.valueOf(var28)) + "m\n");
                  String locus;
                  Iterator var29;
                  if(var28 == 0) {
                     var29 = this.Working_Loci.iterator();

                     while(var29.hasNext()) {
                        locus = (String)var29.next();
                        var26.write(locus + " = " + ((HashMap)this.n_locus_value.get(Integer.valueOf(var28))).get(locus) + "\n");
                     }

                     var26.write("Likelihood: " + this.n_ans.get(Integer.valueOf(var28)) + "\n");
                     var26.write("Probability: " + this.n_probdist.get(Integer.valueOf(var28)) + "\n\n");
                  } else {
                     var29 = this.Working_Loci.iterator();

                     while(true) {
                        while(var29.hasNext()) {
                           locus = (String)var29.next();
                           ArrayList set_of_alleles;
                           int index;
                           if(!locus.contentEquals("AMEL")) {
                              set_of_alleles = new ArrayList();
                              ArrayList var30 = new ArrayList();
                              Iterator var31 = ((Set)((HashMap)this.n_loci_max_alleles.get(Integer.valueOf(var28))).get(locus)).iterator();

                              while(var31.hasNext()) {
                                 index = ((Integer)var31.next()).intValue();
                                 set_of_alleles.add(Integer.valueOf(index));
                              }

                              var31 = set_of_alleles.iterator();

                              while(var31.hasNext()) {
                                 index = ((Integer)var31.next()).intValue();
                                 if(!((ArrayList)Loci_Alleles.get(locus)).contains(Integer.valueOf(index))) {
                                    int index1 = set_of_alleles.indexOf(Integer.valueOf(index));
                                    set_of_alleles.set(index1, Integer.valueOf(0));
                                 }
                              }

                              var31 = set_of_alleles.iterator();

                              while(var31.hasNext()) {
                                 index = ((Integer)var31.next()).intValue();
                                 var30.add(Double.valueOf((double)index / 10.0D));
                              }

                              var26.write(locus + " = " + ((HashMap)this.n_locus_value.get(Integer.valueOf(var28))).get(locus) + " : " + var30 + "\n");
                           } else {
                              set_of_alleles = new ArrayList();
                              set_of_alleles.addAll((Collection)((HashMap)this.n_amel_max_alleles.get(Integer.valueOf(var28))).get(locus));
                              Iterator i$1 = set_of_alleles.iterator();

                              while(i$1.hasNext()) {
                                 String allele = (String)i$1.next();
                                 if(!((ArrayList)AMEL_Alleles.get(locus)).contains(allele)) {
                                    index = set_of_alleles.indexOf(allele);
                                    set_of_alleles.set(index, String.valueOf(0));
                                 }
                              }

                              var26.write(locus + " = " + ((HashMap)this.n_locus_value.get(Integer.valueOf(var28))).get(locus) + " : " + set_of_alleles + "\n");
                           }
                        }

                        var26.write("Likelihood: " + this.n_ans.get(Integer.valueOf(var28)) + "\n");
                        var26.write("Probability: " + this.n_probdist.get(Integer.valueOf(var28)) + "\n\n");
                        break;
                     }
                  }
               }
            }

            var26.close();
         } catch (Throwable var23) {
            var6 = var23;
            throw var23;
         } finally {
            if(var26 != null) {
               if(var6 != null) {
                  try {
                     var26.close();
                  } catch (Throwable var22) {
                     var6.addSuppressed(var22);
                  }
               } else {
                  var26.close();
               }
            }

         }
      } catch (Exception var25) {
         System.err.println(var25.getMessage());
      }

   }

   public void printNocitResults() {
      BigDecimal check1 = new BigDecimal(0);
      BigDecimal total = new BigDecimal(0);

      int e;
      for(e = 0; e <= this.noc; ++e) {
         total = total.add((BigDecimal)this.n_ans.get(Integer.valueOf(e)));
      }

      if(total.compareTo(check1) == 0) {
         for(e = 0; e <= this.noc; ++e) {
            this.n_probdist.put(Integer.valueOf(e), new BigDecimal(0));
         }
      } else {
         for(e = 0; e <= this.noc; ++e) {
            BigDecimal output_file = ((BigDecimal)this.n_ans.get(Integer.valueOf(e))).divide(total, MathContext.DECIMAL128);
            if(output_file.compareTo(check1) == 0) {
               this.n_probdist.put(Integer.valueOf(e), new BigDecimal(0));
            } else {
               this.n_probdist.put(Integer.valueOf(e), output_file);
            }
         }
      }

      try {
         FileWriter var26 = new FileWriter(this.Output_File_Name);
         BufferedWriter var25 = new BufferedWriter(var26);
         Throwable var5 = null;

         try {
            var25.write("Calibration File: " + this.calibrationPath + "\n");
            var25.write("Frequency File: " + this.Freq_File + "\n");
            var25.write("Sample File: " + this.Sample_File_Name + "\n");
            var25.write("Sample DNA input: " + DNA_Mass + "ng\n");
            var25.write("\n");
            Iterator x2 = this.Loci_not_in_calib_samples.iterator();

            String i$;
            while(x2.hasNext()) {
               i$ = (String)x2.next();
               var25.write("Locus " + i$ + " was not included in the calculation because it is not present in the calibration samples\n");
            }

            x2 = this.Loci_not_in_freq_table.iterator();

            while(x2.hasNext()) {
               i$ = (String)x2.next();
               var25.write("Locus " + i$ + " was not included in the calculation because it is not present in the frequency table\n");
            }

            x2 = this.Loci_no_true.iterator();

            while(x2.hasNext()) {
               i$ = (String)x2.next();
               var25.write("Locus " + i$ + " was not included in the calculation because it had too few data points for allelic peaks\n");
            }

            x2 = this.Loci_no_noise.iterator();

            while(x2.hasNext()) {
               i$ = (String)x2.next();
               var25.write("Locus " + i$ + " was not included in the calculation because it had too few data points for noise peaks\n");
            }

            x2 = Sample_Loci.iterator();

            while(x2.hasNext()) {
               i$ = (String)x2.next();
               if(!this.Loci_RStutter.contains(i$) && !i$.contentEquals("AMEL")) {
                  var25.write("Locus " + i$ + " had too few data points for reverse stutter \n");
               }
            }

            x2 = Sample_Loci.iterator();

            while(x2.hasNext()) {
               i$ = (String)x2.next();
               if(!this.Loci_FStutter.contains(i$) && !i$.contentEquals("AMEL")) {
                  var25.write("Locus " + i$ + " had too few data points for forward stutter \n");
               }
            }

            var25.write("\n");

            for(int var27 = 0; var27 <= this.noc; ++var27) {
               if(this.n_ans.get(Integer.valueOf(var27)) != null) {
                  var25.write("Number of contributors: " + var27 + "\n");
                  var25.write("Time taken: " + this.time_taken.get(Integer.valueOf(var27)) + "m\n");
                  String locus;
                  Iterator var28;
                  if(var27 == 0) {
                     var28 = this.Working_Loci.iterator();

                     while(var28.hasNext()) {
                        locus = (String)var28.next();
                        var25.write(locus + " = " + ((HashMap)this.n_locus_value.get(Integer.valueOf(var27))).get(locus) + "\n");
                     }

                     var25.write("Likelihood: " + this.n_ans.get(Integer.valueOf(var27)) + "\n");
                     var25.write("Probability: " + this.n_probdist.get(Integer.valueOf(var27)) + "\n\n");
                  } else {
                     var28 = this.Working_Loci.iterator();

                     while(true) {
                        while(var28.hasNext()) {
                           locus = (String)var28.next();
                           ArrayList set_of_alleles;
                           int index;
                           if(!locus.contentEquals("AMEL")) {
                              set_of_alleles = new ArrayList();
                              ArrayList var29 = new ArrayList();
                              Iterator var30 = ((Set)((HashMap)this.n_loci_max_alleles.get(Integer.valueOf(var27))).get(locus)).iterator();

                              while(var30.hasNext()) {
                                 index = ((Integer)var30.next()).intValue();
                                 set_of_alleles.add(Integer.valueOf(index));
                              }

                              var30 = set_of_alleles.iterator();

                              while(var30.hasNext()) {
                                 index = ((Integer)var30.next()).intValue();
                                 if(!((ArrayList)Loci_Alleles.get(locus)).contains(Integer.valueOf(index))) {
                                    int index1 = set_of_alleles.indexOf(Integer.valueOf(index));
                                    set_of_alleles.set(index1, Integer.valueOf(0));
                                 }
                              }

                              var30 = set_of_alleles.iterator();

                              while(var30.hasNext()) {
                                 index = ((Integer)var30.next()).intValue();
                                 var29.add(Double.valueOf((double)index / 10.0D));
                              }

                              var25.write(locus + " = " + ((HashMap)this.n_locus_value.get(Integer.valueOf(var27))).get(locus) + " : " + var29 + "\n");
                           } else {
                              set_of_alleles = new ArrayList();
                              set_of_alleles.addAll((Collection)((HashMap)this.n_amel_max_alleles.get(Integer.valueOf(var27))).get(locus));
                              Iterator i$1 = set_of_alleles.iterator();

                              while(i$1.hasNext()) {
                                 String allele = (String)i$1.next();
                                 if(!((ArrayList)AMEL_Alleles.get(locus)).contains(allele)) {
                                    index = set_of_alleles.indexOf(allele);
                                    set_of_alleles.set(index, String.valueOf(0));
                                 }
                              }

                              var25.write(locus + " = " + ((HashMap)this.n_locus_value.get(Integer.valueOf(var27))).get(locus) + " : " + set_of_alleles + "\n");
                           }
                        }

                        var25.write("Likelihood: " + this.n_ans.get(Integer.valueOf(var27)) + "\n");
                        var25.write("Probability: " + this.n_probdist.get(Integer.valueOf(var27)) + "\n\n");
                        break;
                     }
                  }
               }
            }

            var25.close();
         } catch (Throwable var22) {
            var5 = var22;
            throw var22;
         } finally {
            if(var25 != null) {
               if(var5 != null) {
                  try {
                     var25.close();
                  } catch (Throwable var21) {
                     var5.addSuppressed(var21);
                  }
               } else {
                  var25.close();
               }
            }

         }
      } catch (Exception var24) {
         System.err.println(var24.getMessage());
      }

   }

   public void runNocit(File sampleFile, File frequencyFile, File outputFile, File calibrationFile, int maxNumberContributors, double dnaMass, JProgressBar progBar) throws InterruptedException, ExecutionException, FunctionEvaluationException, IOException, OptimizationException {
      double epsilon = 0.001D;
      int[] loci_increments = new int[]{100000, 10000000, 60000000, 120000000, 200000000};
      final int[] amel_increments = { 50000, 1000000, 6000000, 60000000, 120000000 };
      LevenbergMarquardtOptimizer lmo = new LevenbergMarquardtOptimizer();
      lmo.setMaxIterations(1000000);
      int step_size = 1000000;
      BigDecimal check = new BigDecimal(0);
      int num_processors = Runtime.getRuntime().availableProcessors();
      int corePoolSize = num_processors;
      int maxPoolSize = num_processors;
      long keepAliveTime = 5L;
      this.n_probdist = new HashMap();
      this.n_ans = new HashMap();
      this.n_locus_value = new HashMap();
      this.n_loci_max_alleles = new HashMap();
      this.n_amel_max_alleles = new HashMap();
      this.time_taken = new HashMap();
      this.Loci_no_true = new HashSet();
      this.Loci_no_noise = new HashSet();
      this.Loci_RStutter = new HashSet();
      this.Loci_FStutter = new HashSet();
      this.Loci_not_in_freq_table = new HashSet();
      this.Loci_not_in_calib_samples = new HashSet();
      this.Working_Loci = new HashSet();
      String male = "M";
      String female = "F";
      Sexes = new String[]{"M", "F"};
      String[] male_genotype = new String[]{"X", "Y"};
      String[] female_genotype = new String[]{"X", "X"};
      Sexes_Genotype = new HashMap();
      Sexes_Genotype.put(male, male_genotype);
      Sexes_Genotype.put(female, female_genotype);
      this.Freq_File = frequencyFile.getAbsolutePath();
      this.calibrationPath = calibrationFile.getAbsolutePath();
      this.Sample_File_Name = sampleFile.getAbsolutePath();
      this.Output_File_Name = outputFile.getAbsolutePath();
      DNA_Mass = dnaMass;
      this.noc = maxNumberContributors;
      NocitAlgorithm.FreqTable FreqTableObject = new NocitAlgorithm.FreqTable(this.Freq_File);
      LinkedHashMap Freq_Table = FreqTableObject.getFreqTable();
      ArrayList Freq_Loci = FreqTableObject.getFreqLociList();
      Intervals_Map = intervals(Freq_Table, Freq_Loci);
      Iterator pSFObject = Freq_Loci.iterator();

      while(pSFObject.hasNext()) {
         String cdcobject = (String)pSFObject.next();
         Set Calib_DO = ((LinkedHashMap)Intervals_Map.get(cdcobject)).keySet();
         ArrayList Calib_rstutDO = new ArrayList();
         Calib_rstutDO.addAll(Calib_DO);
         Integer Calib_fstutDO = (Integer)Calib_rstutDO.get(Calib_rstutDO.size() - 1);
         Final_Values.put(cdcobject, Calib_fstutDO);
      }

      NocitAlgorithm.parseSampleFile var86 = new NocitAlgorithm.parseSampleFile(this.Sample_File_Name);
      Loci_Peaks = var86.getLociPeaks();
      Loci_Alleles = var86.getLociAlleles();
      AMEL_Peaks = var86.getAMELPeaks();
      AMEL_Alleles = var86.getAMELAlleles();
      Sample_Loci = var86.getSampleLoci();
      NocitAlgorithm.CalibDataCollection var87 = new NocitAlgorithm.CalibDataCollection(this.calibrationPath, Freq_Table);
      HashMap var88 = var87.getDOCalibData();
      HashMap var89 = var87.getrstutDOCalibData();
      HashMap var90 = var87.getfstutDOCalibData();
      HashMap Calib_True = var87.getTrueCalibData();
      HashMap Calib_Noise = var87.getNoiseCalibData();
      HashMap Calib_rstut = var87.getrstutterCalibData();
      HashMap Calib_fstut = var87.getfstutterCalibData();
      HashSet Calib_Loci = var87.getLoci();
      HashSet Calib_Masses = var87.getMasses();
      ArrayList calib_masses = new ArrayList(Calib_Masses);
      Collections.sort(calib_masses);
      HashSet Loci_cannot_work = new HashSet();
      True_Mean_Slope = new HashMap();
      True_Stddev_Slope = new HashMap();
      Noise_Mean_Slope = new HashMap();
      Noise_Stddev_Slope = new HashMap();
      R_Stutter_Mean = new HashMap();
      R_Stutter_Stddev = new HashMap();
      F_Stutter_Mean = new HashMap();
      F_Stutter_Stddev = new HashMap();
      Locus_DO = new HashMap();
      R_Stut_DO = new HashMap();
      F_Stut_DO = new HashMap();
      Iterator n = Calib_Loci.iterator();

      String initial_time;
      while(n.hasNext()) {
         initial_time = (String)n.next();
         ArrayList total_True_Masses = new ArrayList();
         ArrayList Max_Val = new ArrayList();
         ArrayList locus_vals = new ArrayList();
         ArrayList locus_zerovals = new ArrayList();
         ArrayList n_val = new ArrayList();
         ArrayList final_time = new ArrayList();
         ArrayList locus = new ArrayList();
         CurveFitter time = new CurveFitter(lmo);
         CurveFitter not_reached = new CurveFitter(lmo);
         CurveFitter time_ans = new CurveFitter(lmo);
         CurveFitter split = new CurveFitter(lmo);
         CurveFitter iterations = new CurveFitter(lmo);
         CurveFitter rstutstddev_fitter = new CurveFitter(lmo);
         CurveFitter present_value = new CurveFitter(lmo);
         CurveFitter fstutstddev_fitter = new CurveFitter(lmo);
         CurveFitter present_transformed_value = new CurveFitter(lmo);
         CurveFitter rstutdo_fitter = new CurveFitter(lmo);
         CurveFitter prev_value = new CurveFitter(lmo);
         double[] truemean_init = new double[]{3000.0D, 30.0D};
         double[] perc_change = new double[]{1000.0D, 50.0D};
         double[] noisemean_init = new double[]{50.0D, 5.0D};
         double[] locus_answer = new double[]{100.0D, 0.0D};
         double[] futures = new double[]{0.0D, -50.0D, 0.0D};
         double[] tpe = new double[]{0.0D, -50.0D, 0.0D};
         double[] i$ = new double[]{0.0D, -50.0D, 0.0D};
         double[] future = new double[]{0.0D, -50.0D, 0.0D};
         double[] clobject = new double[]{0.0D, -100.0D};
         double[] rstutdo_init = new double[]{0.0D, -30.0D};
         double[] fstutdo_init = new double[]{0.0D, -5.0D};
         Iterator truemean_best = calib_masses.iterator();

         while(truemean_best.hasNext()) {
            double truestddev_best = ((Double)truemean_best.next()).doubleValue();
            Double noisestddev_best = (Double)((ArrayList)((HashMap)Calib_True.get(Double.valueOf(truestddev_best))).get(initial_time)).get(0);
            Double rstutmean_best = (Double)((ArrayList)((HashMap)Calib_True.get(Double.valueOf(truestddev_best))).get(initial_time)).get(1);
            if(!rstutmean_best.isNaN() && rstutmean_best.doubleValue() != 0.0D) {
               total_True_Masses.add(Double.valueOf(truestddev_best));
               time.addObservedPoint(truestddev_best, noisestddev_best.doubleValue());
               not_reached.addObservedPoint(truestddev_best, rstutmean_best.doubleValue());
            }

            Double rstutstddev_best = (Double)((ArrayList)((HashMap)Calib_Noise.get(Double.valueOf(truestddev_best))).get(initial_time)).get(0);
            Double fstutmean_best = (Double)((ArrayList)((HashMap)Calib_Noise.get(Double.valueOf(truestddev_best))).get(initial_time)).get(1);
            if(!fstutmean_best.isNaN() && fstutmean_best.doubleValue() != 0.0D) {
               Max_Val.add(Double.valueOf(truestddev_best));
               time_ans.addObservedPoint(truestddev_best, rstutstddev_best.doubleValue());
               split.addObservedPoint(truestddev_best, fstutmean_best.doubleValue());
            }

            Double fstutstddev_best = (Double)((HashMap)var88.get(Double.valueOf(truestddev_best))).get(initial_time);
            if(!fstutstddev_best.isInfinite()) {
               n_val.add(Double.valueOf(truestddev_best));
               present_transformed_value.addObservedPoint(truestddev_best, fstutstddev_best.doubleValue());
            }

            if(!initial_time.contentEquals("AMEL")) {
               Double do_best = (Double)((HashMap)var89.get(Double.valueOf(truestddev_best))).get(initial_time);
               if(!do_best.isInfinite()) {
                  final_time.add(Double.valueOf(truestddev_best));
                  rstutdo_fitter.addObservedPoint(truestddev_best, do_best.doubleValue());
               }

               Double rstutdo_best = (Double)((ArrayList)((HashMap)Calib_rstut.get(Double.valueOf(truestddev_best))).get(initial_time)).get(0);
               Double fstutdo_best = (Double)((ArrayList)((HashMap)Calib_rstut.get(Double.valueOf(truestddev_best))).get(initial_time)).get(1);
               if(!fstutdo_best.isNaN() && fstutdo_best.doubleValue() != 0.0D) {
                  iterations.addObservedPoint(truestddev_best, rstutdo_best.doubleValue());
                  rstutstddev_fitter.addObservedPoint(truestddev_best, fstutdo_best.doubleValue());
                  locus_vals.add(Double.valueOf(truestddev_best));
               }

               Double fstutstddev = (Double)((HashMap)var90.get(Double.valueOf(truestddev_best))).get(initial_time);
               if(!fstutstddev.isInfinite()) {
                  locus.add(Double.valueOf(truestddev_best));
                  prev_value.addObservedPoint(truestddev_best, fstutstddev.doubleValue());
               }

               Double new_fstutstddev_best = (Double)((ArrayList)((HashMap)Calib_fstut.get(Double.valueOf(truestddev_best))).get(initial_time)).get(0);
               Double fstut_stddevval = (Double)((ArrayList)((HashMap)Calib_fstut.get(Double.valueOf(truestddev_best))).get(initial_time)).get(1);
               if(!fstut_stddevval.isNaN() && fstut_stddevval.doubleValue() != 0.0D) {
                  present_value.addObservedPoint(truestddev_best, new_fstutstddev_best.doubleValue());
                  fstutstddev_fitter.addObservedPoint(truestddev_best, fstut_stddevval.doubleValue());
                  locus_zerovals.add(Double.valueOf(truestddev_best));
               }
            }
         }

         double[] var143 = new double[2];
         double[] var144 = new double[2];
         double[] noisemean_best = new double[2];
         double[] var145 = new double[2];
         double[] var146 = new double[3];
         double[] var147 = new double[3];
         double[] var148 = new double[3];
         double[] var149 = new double[3];
         double[] var150 = new double[]{0.0D, 0.0D};
         double[] var151 = new double[]{0.0D, 0.0D};
         double[] var152 = new double[]{0.0D, 0.0D};
         if(total_True_Masses.size() > 1) {
            var143 = time.fit(new NocitAlgorithm.LinFunction(), truemean_init);
            var144 = not_reached.fit(new NocitAlgorithm.LinFunction(), perc_change);
         } else {
            Loci_cannot_work.add(initial_time);
            this.Loci_no_true.add(initial_time);
         }

         if(Max_Val.size() > 1) {
            noisemean_best = time_ans.fit(new NocitAlgorithm.LinFunction(), noisemean_init);
            var145 = split.fit(new NocitAlgorithm.LinFunction(), locus_answer);
         } else {
            Loci_cannot_work.add(initial_time);
            this.Loci_no_noise.add(initial_time);
         }

         if(n_val.size() > 1) {
            var150 = present_transformed_value.fit(new NocitAlgorithm.ExpFunction(), clobject);
         }

         if(!initial_time.contentEquals("AMEL")) {
            if(locus_vals.size() > 2) {
               var146 = iterations.fit(new NocitAlgorithm.TwoExpFunction(), futures);
               var147 = rstutstddev_fitter.fit(new NocitAlgorithm.TwoExpFunction(), tpe);
               this.Loci_RStutter.add(initial_time);
            }

            if(locus_zerovals.size() > 2) {
               var148 = present_value.fit(new NocitAlgorithm.TwoExpFunction(), i$);
               var149 = fstutstddev_fitter.fit(new NocitAlgorithm.TwoExpFunction(), future);
               this.Loci_FStutter.add(initial_time);
            }

            if(final_time.size() > 1) {
               var151 = rstutdo_fitter.fit(new NocitAlgorithm.ExpFunction(), rstutdo_init);
            }

            if(locus.size() > 1) {
               var152 = prev_value.fit(new NocitAlgorithm.ExpFunction(), fstutdo_init);
            }
         }

         if(!Loci_cannot_work.contains(initial_time)) {
            True_Mean_Slope.put(initial_time, var143);
            True_Stddev_Slope.put(initial_time, var144);
            Noise_Mean_Slope.put(initial_time, noisemean_best);
            Noise_Stddev_Slope.put(initial_time, var145);
            Locus_DO.put(initial_time, var150);
         }

         if(!initial_time.contentEquals("AMEL")) {
            double[] var153;
            double[] var154;
            if(this.Loci_RStutter.contains(initial_time)) {
               if(var146[2] >= 0.0D) {
                  R_Stutter_Mean.put(initial_time, var146);
               } else {
                  var153 = iterations.fit(new NocitAlgorithm.PosTwoExpFunction(), futures);
                  if(var153[2] >= 0.0D) {
                     R_Stutter_Mean.put(initial_time, var153);
                  } else {
                     var154 = new double[]{var153[0], var153[1], 0.0D};
                     R_Stutter_Mean.put(initial_time, var154);
                  }
               }

               if(var147[2] >= 0.0D) {
                  R_Stutter_Stddev.put(initial_time, var147);
               } else {
                  var153 = rstutstddev_fitter.fit(new NocitAlgorithm.PosTwoExpFunction(), tpe);
                  if(var153[2] >= 0.0D) {
                     R_Stutter_Stddev.put(initial_time, var153);
                  } else {
                     var154 = new double[]{var153[0], var153[1], 0.0D};
                     R_Stutter_Stddev.put(initial_time, var154);
                  }
               }

               R_Stut_DO.put(initial_time, var151);
            }

            if(this.Loci_FStutter.contains(initial_time)) {
               if(var148[2] >= 0.0D) {
                  F_Stutter_Mean.put(initial_time, var148);
               } else {
                  var153 = present_value.fit(new NocitAlgorithm.PosTwoExpFunction(), i$);
                  if(var153[2] >= 0.0D) {
                     F_Stutter_Mean.put(initial_time, var153);
                  } else {
                     var154 = new double[]{var153[0], var153[1], 0.0D};
                     F_Stutter_Mean.put(initial_time, var154);
                  }
               }

               if(var149[2] >= 0.0D) {
                  F_Stutter_Stddev.put(initial_time, var149);
               } else {
                  var153 = fstutstddev_fitter.fit(new NocitAlgorithm.PosTwoExpFunction(), future);
                  if(var153[2] >= 0.0D) {
                     F_Stutter_Stddev.put(initial_time, var153);
                  } else {
                     var154 = new double[]{var153[0], var153[1], 0.0D};
                     F_Stutter_Stddev.put(initial_time, var154);
                  }
               }

               F_Stut_DO.put(initial_time, var152);
            }
         }
      }

      n = Sample_Loci.iterator();

      while(n.hasNext()) {
         initial_time = (String)n.next();
         if(Calib_Loci.contains(initial_time)) {
            if(!Loci_cannot_work.contains(initial_time)) {
               if(!initial_time.contentEquals("AMEL")) {
                  if(Freq_Loci.contains(initial_time)) {
                     this.Working_Loci.add(initial_time);
                  } else {
                     this.Loci_not_in_freq_table.add(initial_time);
                  }
               } else {
                  this.Working_Loci.add(initial_time);
               }
            }
         } else {
            this.Loci_not_in_calib_samples.add(initial_time);
         }
      }

      for(int var91 = 0; var91 <= this.noc; ++var91) {
         double var92;
         HashMap var99;
         BigDecimal var103;
         double var108;
         double var122;
         double var123;
         if(var91 == 0) {
            var92 = (double)System.currentTimeMillis();
            this.n_locus_value.put(Integer.valueOf(var91), new HashMap());
            Iterator var94 = this.Working_Loci.iterator();

            while(true) {
               while(var94.hasNext()) {
                  String var97 = (String)var94.next();
                  HashMap var105;
                  ConcurrentHashMap var109;
                  Iterator var116;
                  double[] var119;
                  BigDecimal var120;
                  if(!var97.contentEquals("AMEL")) {
                     var99 = new HashMap();
                     var105 = new HashMap();
                     var109 = (ConcurrentHashMap)Loci_Peaks.get(var97);
                     var116 = var109.values().iterator();

                     while(var116.hasNext()) {
                        NocitAlgorithm.Peak var117 = (NocitAlgorithm.Peak)var116.next();
                        var119 = calcSlopeValue(var97, DNA_Mass, Noise_Mean_Slope, Noise_Stddev_Slope);
                        var123 = var119[0];
                        var122 = var119[1] * var119[1];
                        var99.put(Integer.valueOf(var117.getAllele()), Double.valueOf(var123));
                        var105.put(Integer.valueOf(var117.getAllele()), Double.valueOf(var122));
                     }

                     var108 = calcLocusPeakHeightsProb(var109, var99, var105);
                     var120 = new BigDecimal(var108, MathContext.DECIMAL128);
                     ((HashMap)this.n_locus_value.get(Integer.valueOf(var91))).put(var97, var120);
                  } else {
                     var99 = new HashMap();
                     var105 = new HashMap();
                     var109 = (ConcurrentHashMap)AMEL_Peaks.get(var97);
                     var116 = var109.values().iterator();

                     while(var116.hasNext()) {
                        NocitAlgorithm.AMEL_Peak var114 = (NocitAlgorithm.AMEL_Peak)var116.next();
                        var119 = calcSlopeValue(var97, DNA_Mass, Noise_Mean_Slope, Noise_Stddev_Slope);
                        var123 = var119[0];
                        var122 = var119[1] * var119[1];
                        var99.put(var114.getAllele(), Double.valueOf(var123));
                        var105.put(var114.getAllele(), Double.valueOf(var122));
                     }

                     var108 = calcAMELPeakHeightsProb(var109, var99, var105);
                     var120 = new BigDecimal(var108, MathContext.DECIMAL128);
                     ((HashMap)this.n_locus_value.get(Integer.valueOf(var91))).put(var97, var120);
                  }
               }

               BigDecimal var95 = new BigDecimal(1, MathContext.DECIMAL128);

               for(Iterator var98 = this.Working_Loci.iterator(); var98.hasNext(); var95 = var95.multiply(var103, MathContext.DECIMAL128)) {
                  String var102 = (String)var98.next();
                  var103 = (BigDecimal)((HashMap)this.n_locus_value.get(Integer.valueOf(var91))).get(var102);
               }

               if(var95.compareTo(check) == 0) {
                  this.n_ans.put(Integer.valueOf(var91), new BigDecimal(0));
               } else {
                  this.n_ans.put(Integer.valueOf(var91), var95);
               }

               double var100 = (double)System.currentTimeMillis();
               double var110 = (var100 - var92) / 60000.0D;
               var108 = (double)Math.round(var110 * 100.0D) / 100.0D;
               this.time_taken.put(Integer.valueOf(var91), Double.valueOf(var108));
               break;
            }
         } else {
            var92 = (double)System.currentTimeMillis();
            this.n_loci_max_alleles.put(Integer.valueOf(var91), new HashMap());
            this.n_amel_max_alleles.put(Integer.valueOf(var91), new HashMap());
            this.n_locus_value.put(Integer.valueOf(var91), new HashMap());
            HashMap var93 = new HashMap();
            HashMap var96 = new HashMap();
            var99 = new HashMap();
            Iterator var101 = this.Working_Loci.iterator();

            while(var101.hasNext()) {
               String var104 = (String)var101.next();
               var93.put(var104, Double.valueOf(0.0D));
               var96.put(var104, new ArrayList());
               var99.put(var104, new ArrayList());
               var108 = 0.0D;
               boolean var115 = true;
               int var118;
               int var121;
               if(!var104.contentEquals("AMEL")) {
                  var118 = loci_increments[var91 - 1];
                  var121 = var118 / num_processors;
                  var122 = (double)loci_increments[var91 - 1];
               } else {
                  var118 = amel_increments[var91 - 1];
                  var121 = var118 / num_processors;
                  var122 = (double)amel_increments[var91 - 1];
               }

               while(var115) {
                  int var128;
                  ExecutorService var129;
                  ArrayList var130;
                  HashSet var131;
                  ThreadPoolExecutor var132;
                  int var134;
                  List var135;
                  int var136;
                  Iterator var137;
                  int var138;
                  Future var139;
                  Future var140;
                  if(!var104.contentEquals("AMEL")) {
                     NocitAlgorithm.Callable_loci_object var142;
                     if(var91 < 2) {
                        var129 = Executors.newFixedThreadPool(num_processors);
                        var131 = new HashSet();

                        for(var134 = 1; var134 <= num_processors; ++var134) {
                           var131.add(new NocitAlgorithm.Sum_threads_stutterloci(var121, var91, var104, this.Loci_RStutter, this.Loci_FStutter));
                        }

                        var135 = var129.invokeAll(var131);
                        var137 = var135.iterator();

                        while(var137.hasNext()) {
                           var139 = (Future)var137.next();
                           var142 = (NocitAlgorithm.Callable_loci_object)var139.get();
                           var108 += var142.getvalue();
                           if(var142.getmaxvalue() > ((Double)var93.get(var104)).doubleValue()) {
                              var93.put(var104, Double.valueOf(var142.getmaxvalue()));
                              ((HashMap)this.n_loci_max_alleles.get(Integer.valueOf(var91))).put(var104, var142.getmaxalls());
                           }
                        }

                        var129.shutdown();
                     } else {
                        var128 = 0;
                        var130 = new ArrayList();
                        var132 = new ThreadPoolExecutor(corePoolSize, maxPoolSize, keepAliveTime, TimeUnit.MINUTES, new LinkedBlockingQueue());

                        label337:
                        while(true) {
                           if(var128 >= var118) {
                              var137 = var130.iterator();

                              while(true) {
                                 if(!var137.hasNext()) {
                                    break label337;
                                 }

                                 var139 = (Future)var137.next();
                                 var142 = (NocitAlgorithm.Callable_loci_object)var139.get();
                                 var108 += var142.getvalue();
                                 if(var142.getmaxvalue() > ((Double)var93.get(var104)).doubleValue()) {
                                    var93.put(var104, Double.valueOf(var142.getmaxvalue()));
                                    ((HashMap)this.n_loci_max_alleles.get(Integer.valueOf(var91))).put(var104, var142.getmaxalls());
                                 }
                              }
                           }

                           var136 = num_processors - var132.getActiveCount();

                           for(var138 = 1; var138 <= var136; ++var138) {
                              if(var128 < var118) {
                                 var128 += step_size;
                                 var140 = var132.submit(new NocitAlgorithm.Sum_threads_stutterloci(step_size, var91, var104, this.Loci_RStutter, this.Loci_FStutter));
                                 var130.add(var140);
                              }
                           }
                        }
                     }
                  } else {
                     NocitAlgorithm.Callable_amel_object var141;
                     if(var91 <= 3) {
                        var129 = Executors.newFixedThreadPool(num_processors);
                        var131 = new HashSet();

                        for(var134 = 1; var134 <= num_processors; ++var134) {
                           var131.add(new NocitAlgorithm.Sum_threads_amel(var121, var91, var104));
                        }

                        var135 = var129.invokeAll(var131);
                        var137 = var135.iterator();

                        while(var137.hasNext()) {
                           var139 = (Future)var137.next();
                           var141 = (NocitAlgorithm.Callable_amel_object)var139.get();
                           var108 += var141.getvalue();
                           if(var141.getmaxvalue() > ((Double)var93.get(var104)).doubleValue()) {
                              var93.put(var104, Double.valueOf(var141.getmaxvalue()));
                              ((HashMap)this.n_amel_max_alleles.get(Integer.valueOf(var91))).put(var104, var141.getmaxalls());
                           }
                        }

                        var129.shutdown();
                     } else {
                        var128 = 0;
                        var130 = new ArrayList();
                        var132 = new ThreadPoolExecutor(corePoolSize, maxPoolSize, keepAliveTime, TimeUnit.MINUTES, new LinkedBlockingQueue());

                        label361:
                        while(true) {
                           if(var128 >= var118) {
                              var137 = var130.iterator();

                              while(true) {
                                 if(!var137.hasNext()) {
                                    break label361;
                                 }

                                 var139 = (Future)var137.next();
                                 var141 = (NocitAlgorithm.Callable_amel_object)var139.get();
                                 var108 += var141.getvalue();
                                 if(var141.getmaxvalue() > ((Double)var93.get(var104)).doubleValue()) {
                                    var93.put(var104, Double.valueOf(var141.getmaxvalue()));
                                    ((HashMap)this.n_amel_max_alleles.get(Integer.valueOf(var91))).put(var104, var141.getmaxalls());
                                 }
                              }
                           }

                           var136 = num_processors - var132.getActiveCount();

                           for(var138 = 1; var138 <= var136; ++var138) {
                              if(var128 < var118) {
                                 var128 += step_size;
                                 var140 = var132.submit(new NocitAlgorithm.Sum_threads_amel(step_size, var91, var104));
                                 var130.add(var140);
                              }
                           }
                        }
                     }
                  }

                  double var124 = var108 / var122;
                  double var125 = Math.log10(var124);
                  if(var125 != Double.NEGATIVE_INFINITY) {
                     ((ArrayList)var96.get(var104)).add(Double.valueOf(var125));
                     if(((ArrayList)var96.get(var104)).size() > 1) {
                        double var126 = ((Double)((ArrayList)var96.get(var104)).get(((ArrayList)var96.get(var104)).size() - 2)).doubleValue();
                        double var127 = Math.abs((var125 - var126) / var126);
                        if(var127 <= epsilon) {
                           var115 = false;
                           BigDecimal var133 = new BigDecimal(var124, MathContext.DECIMAL128);
                           ((HashMap)this.n_locus_value.get(Integer.valueOf(var91))).put(var104, var133);
                        } else {
                           var122 += (double)var118;
                        }
                     } else {
                        var122 += (double)var118;
                     }
                  } else {
                     ((ArrayList)var99.get(var104)).add(Double.valueOf(var125));
                     if(((ArrayList)var99.get(var104)).size() > 1) {
                        var115 = false;
                        ((HashMap)this.n_locus_value.get(Integer.valueOf(var91))).put(var104, new BigDecimal(0));
                     } else {
                        var122 += (double)var118;
                     }
                  }
               }
            }

            var103 = new BigDecimal(1, MathContext.DECIMAL128);
            Iterator var106 = this.Working_Loci.iterator();

            while(var106.hasNext()) {
               String var113 = (String)var106.next();
               BigDecimal var111 = (BigDecimal)((HashMap)this.n_locus_value.get(Integer.valueOf(var91))).get(var113);
               var103 = var103.multiply(var111, MathContext.DECIMAL128);
               if(var113.contentEquals("AMEL")) {
                  if(((HashMap)this.n_amel_max_alleles.get(Integer.valueOf(var91))).get(var113) == null) {
                     ((HashMap)this.n_amel_max_alleles.get(Integer.valueOf(var91))).put(var113, new HashSet());
                  }
               } else if(((HashMap)this.n_loci_max_alleles.get(Integer.valueOf(var91))).get(var113) == null) {
                  ((HashMap)this.n_loci_max_alleles.get(Integer.valueOf(var91))).put(var113, new HashSet());
               }
            }

            if(var103.compareTo(check) == 0) {
               this.n_ans.put(Integer.valueOf(var91), new BigDecimal(0));
            } else {
               this.n_ans.put(Integer.valueOf(var91), var103);
            }

            double var107 = (double)System.currentTimeMillis();
            double var112 = (var107 - var92) / 60000.0D;
            var123 = (double)Math.round(var112 * 100.0D) / 100.0D;
            this.time_taken.put(Integer.valueOf(var91), Double.valueOf(var123));
         }
      }

   }

   // $FF: synthetic method
   static double[] access$300(String x0, double x1, HashMap x2, HashMap x3) {
      return calcSlopeValue(x0, x1, x2, x3);
   }

   // $FF: synthetic class
   static class SyntheticClass_1 {
   }

   private static class Sum_threads_stutterloci implements Callable {
      private int count_num;
      private int num_of_ppl;
      private String locus_name;
      private double sum;
      private double max_val;
      private Set max_alls = new HashSet();
      private HashSet rstut_loci;
      private HashSet fstut_loci;

      Sum_threads_stutterloci(int count, int numofppl, String locusname, HashSet rstutloci, HashSet fstutloci) {
         this.count_num = count;
         this.num_of_ppl = numofppl;
         this.locus_name = locusname;
         this.sum = 0.0D;
         this.max_val = 0.0D;
         this.rstut_loci = rstutloci;
         this.fstut_loci = fstutloci;
      }

      public NocitAlgorithm.Callable_loci_object call() {
         Integer final_value = (Integer)NocitAlgorithm.Final_Values.get(this.locus_name);
         ArrayList loc_alleles = (ArrayList)NocitAlgorithm.Loci_Alleles.get(this.locus_name);
         HashMap genotypes = new HashMap();
         HashMap weights = new HashMap();
         HashMap means = new HashMap();
         HashMap variances = new HashMap();
         HashSet selected_true_alleles = new HashSet();
         HashSet true_alleles = new HashSet();
         double[] values = NocitAlgorithm.access$300(this.locus_name, NocitAlgorithm.DNA_Mass, NocitAlgorithm.Noise_Mean_Slope, NocitAlgorithm.Noise_Stddev_Slope);
         double noise_mu = values[0];
         double noise_sigmasquared = values[1] * values[1];

         for(int clobject = 1; clobject <= this.count_num; ++clobject) {
            ConcurrentHashMap allPeaks = new ConcurrentHashMap((Map)NocitAlgorithm.Loci_Peaks.get(this.locus_name));
            double[] g;
            if(this.num_of_ppl == 1) {
               g = new double[]{1.0D};
            } else {
               g = NocitAlgorithm.generateRandomTotal(this.num_of_ppl);
            }

            NocitAlgorithm.Peak p;
            int allele_name;
            for(int i$ = 1; i$ <= this.num_of_ppl; ++i$) {
               genotypes.put(Integer.valueOf(i$), new ArrayList());

               for(allele_name = 1; allele_name <= 2; ++allele_name) {
                  int r = (int)Math.round(ThreadLocalRandom.current().nextDouble(1.0D) * (double)final_value.intValue());
                  int allele = ((Integer)((LinkedHashMap)NocitAlgorithm.Intervals_Map.get(this.locus_name)).get(Integer.valueOf(r))).intValue();
                  ((ArrayList)genotypes.get(Integer.valueOf(i$))).add(Integer.valueOf(allele));
                  selected_true_alleles.add(Integer.valueOf(allele));
                  NocitAlgorithm.Peak parentWeight = (NocitAlgorithm.Peak)allPeaks.get(Integer.valueOf(allele));
                  if(parentWeight == null) {
                     p = new NocitAlgorithm.Peak(allele, 0);
                     allPeaks.put(Integer.valueOf(allele), p);
                  }

                  double cont_mass = g[i$ - 1] * NocitAlgorithm.DNA_Mass;
                  if(ThreadLocalRandom.current().nextDouble(1.0D) > NocitAlgorithm.calcExpValue(this.locus_name, NocitAlgorithm.Locus_DO, cont_mass)) {
                     true_alleles.add(Integer.valueOf(allele));
                     if(weights.containsKey(Integer.valueOf(allele))) {
                        weights.put(Integer.valueOf(allele), Double.valueOf(((Double)weights.get(Integer.valueOf(allele))).doubleValue() + cont_mass));
                     } else {
                        weights.put(Integer.valueOf(allele), Double.valueOf(cont_mass));
                     }
                  }
               }
            }

            Iterator var45 = true_alleles.iterator();

            while(var45.hasNext()) {
               allele_name = ((Integer)var45.next()).intValue();
               int height = ((NocitAlgorithm.Peak)allPeaks.get(Integer.valueOf(allele_name))).getHeight();
               values = NocitAlgorithm.access$300(this.locus_name, ((Double)weights.get(Integer.valueOf(allele_name))).doubleValue(), NocitAlgorithm.True_Mean_Slope, NocitAlgorithm.True_Stddev_Slope);
               double true_mu = values[0];
               double true_sigmasquared = values[1] * values[1];
               means.put(Integer.valueOf(allele_name), Double.valueOf(true_mu));
               variances.put(Integer.valueOf(allele_name), Double.valueOf(true_sigmasquared));
               if(loc_alleles.contains(Integer.valueOf(allele_name))) {
                  double var47 = ((Double)weights.get(Integer.valueOf(allele_name))).doubleValue();
                  NocitAlgorithm.Peak fowStutterPeak;
                  if(this.rstut_loci.contains(this.locus_name) && ThreadLocalRandom.current().nextDouble(1.0D) > NocitAlgorithm.calcExpValue(this.locus_name, NocitAlgorithm.R_Stut_DO, var47)) {
                     double rstut_mu = NocitAlgorithm.calcTwoExpValue(this.locus_name, NocitAlgorithm.R_Stutter_Mean, var47) * (double)height;
                     double rstut_sigma = NocitAlgorithm.calcTwoExpValue(this.locus_name, NocitAlgorithm.R_Stutter_Stddev, var47) * (double)height;
                     int revStutterAllele = allele_name - 10;
                     fowStutterPeak = (NocitAlgorithm.Peak)allPeaks.get(Integer.valueOf(revStutterAllele));
                     if(fowStutterPeak == null) {
                        p = new NocitAlgorithm.Peak(revStutterAllele, 0);
                        allPeaks.put(Integer.valueOf(revStutterAllele), p);
                     }

                     if(means.containsKey(Integer.valueOf(revStutterAllele))) {
                        means.put(Integer.valueOf(revStutterAllele), Double.valueOf(((Double)means.get(Integer.valueOf(revStutterAllele))).doubleValue() + rstut_mu));
                        variances.put(Integer.valueOf(revStutterAllele), Double.valueOf(((Double)variances.get(Integer.valueOf(revStutterAllele))).doubleValue() + rstut_sigma * rstut_sigma));
                     } else {
                        means.put(Integer.valueOf(revStutterAllele), Double.valueOf(rstut_mu));
                        variances.put(Integer.valueOf(revStutterAllele), Double.valueOf(rstut_sigma * rstut_sigma));
                     }
                  }

                  if(this.fstut_loci.contains(this.locus_name) && ThreadLocalRandom.current().nextDouble(1.0D) > NocitAlgorithm.calcExpValue(this.locus_name, NocitAlgorithm.F_Stut_DO, var47)) {
                     double fstut_mu = NocitAlgorithm.calcTwoExpValue(this.locus_name, NocitAlgorithm.F_Stutter_Mean, var47) * (double)height;
                     double fstut_sigma = NocitAlgorithm.calcTwoExpValue(this.locus_name, NocitAlgorithm.F_Stutter_Stddev, var47) * (double)height;
                     int fowStutterAllele = allele_name + 10;
                     fowStutterPeak = (NocitAlgorithm.Peak)allPeaks.get(Integer.valueOf(fowStutterAllele));
                     if(fowStutterPeak == null) {
                        p = new NocitAlgorithm.Peak(fowStutterAllele, 0);
                        allPeaks.put(Integer.valueOf(fowStutterAllele), p);
                     }

                     if(means.containsKey(Integer.valueOf(fowStutterAllele))) {
                        means.put(Integer.valueOf(fowStutterAllele), Double.valueOf(((Double)means.get(Integer.valueOf(fowStutterAllele))).doubleValue() + fstut_mu));
                        variances.put(Integer.valueOf(fowStutterAllele), Double.valueOf(((Double)variances.get(Integer.valueOf(fowStutterAllele))).doubleValue() + fstut_sigma * fstut_sigma));
                     } else {
                        means.put(Integer.valueOf(fowStutterAllele), Double.valueOf(fstut_mu));
                        variances.put(Integer.valueOf(fowStutterAllele), Double.valueOf(fstut_sigma * fstut_sigma));
                     }
                  }
               }
            }

            var45 = allPeaks.values().iterator();

            while(var45.hasNext()) {
               NocitAlgorithm.Peak var46 = (NocitAlgorithm.Peak)var45.next();
               if(!means.containsKey(Integer.valueOf(var46.getAllele()))) {
                  means.put(Integer.valueOf(var46.getAllele()), Double.valueOf(noise_mu));
                  variances.put(Integer.valueOf(var46.getAllele()), Double.valueOf(noise_sigmasquared));
               }
            }

            double val_y = NocitAlgorithm.calcLocusPeakHeightsProb(allPeaks, means, variances);
            this.sum += val_y;
            if(val_y > this.max_val) {
               this.max_val = val_y;
               this.max_alls.clear();
               var45 = selected_true_alleles.iterator();

               while(var45.hasNext()) {
                  allele_name = ((Integer)var45.next()).intValue();
                  this.max_alls.add(Integer.valueOf(allele_name));
               }
            }

            genotypes.clear();
            weights.clear();
            means.clear();
            variances.clear();
            true_alleles.clear();
            selected_true_alleles.clear();
            allPeaks.clear();
         }

         NocitAlgorithm.Callable_loci_object var44 = new NocitAlgorithm.Callable_loci_object(this.sum, this.max_val, this.max_alls);
         return var44;
      }
   }

   private static class Sum_threads_amel implements Callable {
      private int count_num;
      private int num_of_ppl;
      private String locus_name;
      private double max_val;
      private Set max_alls = new HashSet();
      private double sum;

      Sum_threads_amel(int count, int numofppl, String locusname) {
         this.count_num = count;
         this.num_of_ppl = numofppl;
         this.locus_name = locusname;
         this.max_val = 0.0D;
         this.sum = 0.0D;
      }

      public NocitAlgorithm.Callable_amel_object call() {
         HashMap genotypes = new HashMap();
         HashMap weights = new HashMap();
         HashMap means = new HashMap();
         HashMap variances = new HashMap();
         HashSet true_alleles = new HashSet();
         HashSet selected_true_alleles = new HashSet();
         double[] values = NocitAlgorithm.access$300(this.locus_name, NocitAlgorithm.DNA_Mass, NocitAlgorithm.Noise_Mean_Slope, NocitAlgorithm.Noise_Stddev_Slope);
         double noise_mu = values[0];
         double noise_sigmasquared = values[1] * values[1];

         for(int clobject = 1; clobject <= this.count_num; ++clobject) {
            ConcurrentHashMap allPeaks = new ConcurrentHashMap((Map)NocitAlgorithm.AMEL_Peaks.get(this.locus_name));
            double[] g;
            if(this.num_of_ppl == 1) {
               g = new double[]{1.0D};
            } else {
               g = NocitAlgorithm.generateRandomTotal(this.num_of_ppl);
            }

            String allele_name;
            for(int i$ = 1; i$ <= this.num_of_ppl; ++i$) {
               genotypes.put(Integer.valueOf(i$), new ArrayList(2));
               allele_name = NocitAlgorithm.Sexes[ThreadLocalRandom.current().nextInt(2)];
               ((ArrayList)genotypes.get(Integer.valueOf(i$))).add(((String[])NocitAlgorithm.Sexes_Genotype.get(allele_name))[0]);
               ((ArrayList)genotypes.get(Integer.valueOf(i$))).add(((String[])NocitAlgorithm.Sexes_Genotype.get(allele_name))[1]);
               Iterator i$1 = ((ArrayList)genotypes.get(Integer.valueOf(i$))).iterator();

               while(i$1.hasNext()) {
                  String allele = (String)i$1.next();
                  NocitAlgorithm.AMEL_Peak peakobj = (NocitAlgorithm.AMEL_Peak)allPeaks.get(allele);
                  if(peakobj == null) {
                     NocitAlgorithm.AMEL_Peak p = new NocitAlgorithm.AMEL_Peak(allele, 0);
                     allPeaks.put(allele, p);
                  }

                  selected_true_alleles.add(allele);
                  double cont_mass = g[i$ - 1] * NocitAlgorithm.DNA_Mass;
                  if(ThreadLocalRandom.current().nextDouble(1.0D) > NocitAlgorithm.calcExpValue(this.locus_name, NocitAlgorithm.Locus_DO, cont_mass)) {
                     true_alleles.add(allele);
                     if(weights.containsKey(allele)) {
                        weights.put(allele, Double.valueOf(((Double)weights.get(allele)).doubleValue() + cont_mass));
                     } else {
                        weights.put(allele, Double.valueOf(cont_mass));
                     }
                  }
               }
            }

            Iterator var30 = true_alleles.iterator();

            while(var30.hasNext()) {
               allele_name = (String)var30.next();
               values = NocitAlgorithm.access$300(this.locus_name, ((Double)weights.get(allele_name)).doubleValue(), NocitAlgorithm.True_Mean_Slope, NocitAlgorithm.True_Stddev_Slope);
               double true_mu = values[0];
               double true_sigmasquared = values[1] * values[1];
               means.put(allele_name, Double.valueOf(true_mu));
               variances.put(allele_name, Double.valueOf(true_sigmasquared));
            }

            var30 = allPeaks.values().iterator();

            while(var30.hasNext()) {
               NocitAlgorithm.AMEL_Peak var31 = (NocitAlgorithm.AMEL_Peak)var30.next();
               if(!means.containsKey(var31.getAllele())) {
                  means.put(var31.getAllele(), Double.valueOf(noise_mu));
                  variances.put(var31.getAllele(), Double.valueOf(noise_sigmasquared));
               }
            }

            double val_y = NocitAlgorithm.calcAMELPeakHeightsProb(allPeaks, means, variances);
            this.sum += val_y;
            if(val_y > this.max_val) {
               this.max_val = val_y;
               this.max_alls.clear();
               var30 = selected_true_alleles.iterator();

               while(var30.hasNext()) {
                  allele_name = (String)var30.next();
                  this.max_alls.add(allele_name);
               }
            }

            genotypes.clear();
            weights.clear();
            means.clear();
            variances.clear();
            true_alleles.clear();
            allPeaks.clear();
            selected_true_alleles.clear();
         }

         NocitAlgorithm.Callable_amel_object var29 = new NocitAlgorithm.Callable_amel_object(this.sum, this.max_val, this.max_alls);
         return var29;
      }
   }

   private static class parseSampleFile {
      private String sample_file_path;
      private ConcurrentHashMap loci_peaks;
      private ConcurrentHashMap AMEL_peaks;
      private HashMap loci_peakalleles;
      private HashMap AMEL_peakalleles;
      private ArrayList sample_loci;

      parseSampleFile(String samplefilepath) {
         this.sample_file_path = samplefilepath;
         this.loci_peaks = new ConcurrentHashMap();
         this.AMEL_peaks = new ConcurrentHashMap();
         this.loci_peakalleles = new HashMap();
         this.AMEL_peakalleles = new HashMap();
         this.sample_loci = new ArrayList();

         try {
            NocitAlgorithm.OpenAndReadFile e = new NocitAlgorithm.OpenAndReadFile(this.sample_file_path);
            String[] samplefilelines = e.FileContents();

            label150:
            for(int linecount = 1; linecount < samplefilelines.length; ++linecount) {
               String line = samplefilelines[linecount];
               String[] line_parts = line.split(",");
               ArrayList line_refresh = new ArrayList();
               String[] locus = line_parts;
               int peakalleles_deleted = line_parts.length;

               for(int peaks_tobedeleted = 0; peaks_tobedeleted < peakalleles_deleted; ++peaks_tobedeleted) {
                  String peak_array = locus[peaks_tobedeleted];
                  if(peak_array != null) {
                     line_refresh.add(peak_array);
                  }
               }

               String var27 = (String)line_refresh.get(1);
               this.sample_loci.add(var27);
               ArrayList peak_alleles;
               ArrayList peak_alleles_refresh;
               int i;
               int count;
               Iterator i$1;
               int index1;
               int index2;
               int h1;
               int h2;
               ArrayList var28;
               ArrayList var29;
               ArrayList var30;
               Iterator var31;
               if(!var27.equals("AMEL")) {
                  this.loci_peaks.put(var27, new ConcurrentHashMap());
                  var28 = new ArrayList();
                  var29 = new ArrayList();
                  var30 = new ArrayList();
                  peak_alleles = new ArrayList();
                  peak_alleles_refresh = new ArrayList();

                  int var36;
                  for(i = 3; i < line_refresh.size(); i += 3) {
                     if(!((String)line_refresh.get(i)).equals("OL")) {
                        double var32 = Double.parseDouble((String)line_refresh.get(i));
                        var36 = (int)Math.round(var32 * 10.0D);
                        count = Integer.parseInt((String)line_refresh.get(i + 2));
                        NocitAlgorithm.Peak var37 = new NocitAlgorithm.Peak(var36, count);
                        var30.add(var37);
                        peak_alleles.add(Integer.valueOf(var37.getAllele()));
                     }
                  }

                  var31 = var30.iterator();

                  while(true) {
                     NocitAlgorithm.Peak var34;
                     do {
                        if(!var31.hasNext()) {
                           var31 = var30.iterator();

                           while(var31.hasNext()) {
                              var34 = (NocitAlgorithm.Peak)var31.next();
                              if(!var29.contains(var34)) {
                                 peak_alleles_refresh.add(Integer.valueOf(var34.getAllele()));
                                 ((ConcurrentHashMap)this.loci_peaks.get(var27)).put(Integer.valueOf(var34.getAllele()), var34);
                              }
                           }

                           this.loci_peakalleles.put(var27, peak_alleles_refresh);
                           continue label150;
                        }

                        var34 = (NocitAlgorithm.Peak)var31.next();
                        var36 = var34.getAllele();
                        count = Collections.frequency(peak_alleles, Integer.valueOf(var36));
                     } while(count <= 1);

                     i$1 = var30.iterator();

                     while(i$1.hasNext()) {
                        NocitAlgorithm.Peak var38 = (NocitAlgorithm.Peak)i$1.next();
                        int var39 = var38.getAllele();
                        if(var36 == var39) {
                           index1 = var30.indexOf(var34);
                           index2 = var30.indexOf(var38);
                           if(index1 != index2) {
                              h1 = var34.getHeight();
                              h2 = var38.getHeight();
                              if(h1 >= h2) {
                                 if(!var28.contains(Integer.valueOf(var39))) {
                                    var29.add(var38);
                                    var28.add(Integer.valueOf(var39));
                                 }
                              } else if(!var28.contains(Integer.valueOf(var36))) {
                                 var29.add(var34);
                                 var28.add(Integer.valueOf(var36));
                              }
                           }
                        }
                     }
                  }
               } else {
                  this.AMEL_peaks.put(var27, new ConcurrentHashMap());
                  var28 = new ArrayList();
                  var29 = new ArrayList();
                  var30 = new ArrayList();
                  peak_alleles = new ArrayList();
                  peak_alleles_refresh = new ArrayList();

                  for(i = 3; i < line_refresh.size(); i += 3) {
                     if(!((String)line_refresh.get(i)).equals("OL")) {
                        String i$ = (String)line_refresh.get(i);
                        int peakobj = Integer.parseInt((String)line_refresh.get(i + 2));
                        NocitAlgorithm.AMEL_Peak allele1 = new NocitAlgorithm.AMEL_Peak(i$, peakobj);
                        var30.add(allele1);
                        peak_alleles.add(allele1.getAllele());
                     }
                  }

                  var31 = var30.iterator();

                  while(true) {
                     NocitAlgorithm.AMEL_Peak var33;
                     String var35;
                     do {
                        if(!var31.hasNext()) {
                           var31 = var30.iterator();

                           while(var31.hasNext()) {
                              var33 = (NocitAlgorithm.AMEL_Peak)var31.next();
                              if(!var29.contains(var33)) {
                                 peak_alleles_refresh.add(var33.getAllele());
                                 ((ConcurrentHashMap)this.AMEL_peaks.get(var27)).put(var33.getAllele(), var33);
                              }
                           }

                           this.AMEL_peakalleles.put(var27, peak_alleles_refresh);
                           continue label150;
                        }

                        var33 = (NocitAlgorithm.AMEL_Peak)var31.next();
                        var35 = var33.getAllele();
                        count = Collections.frequency(peak_alleles, var35);
                     } while(count <= 1);

                     i$1 = var30.iterator();

                     while(i$1.hasNext()) {
                        NocitAlgorithm.AMEL_Peak peakobj2 = (NocitAlgorithm.AMEL_Peak)i$1.next();
                        String allele2 = peakobj2.getAllele();
                        if(var35.contentEquals(allele2)) {
                           index1 = var30.indexOf(var33);
                           index2 = var30.indexOf(peakobj2);
                           if(index1 != index2) {
                              h1 = var33.getHeight();
                              h2 = peakobj2.getHeight();
                              if(h1 >= h2) {
                                 if(!var28.contains(allele2)) {
                                    var29.add(peakobj2);
                                    var28.add(allele2);
                                 }
                              } else if(!var28.contains(var35)) {
                                 var29.add(var33);
                                 var28.add(var35);
                              }
                           }
                        }
                     }
                  }
               }
            }
         } catch (IOException var26) {
            System.out.println(var26.getMessage());
         }

      }

      public HashMap getAMELAlleles() {
         return this.AMEL_peakalleles;
      }

      public ConcurrentHashMap getAMELPeaks() {
         return this.AMEL_peaks;
      }

      public HashMap getLociAlleles() {
         return this.loci_peakalleles;
      }

      public ConcurrentHashMap getLociPeaks() {
         return this.loci_peaks;
      }

      public ArrayList getSampleLoci() {
         return this.sample_loci;
      }
   }

   private static class OpenAndReadFile {
      private String path;

      OpenAndReadFile(String file_path) {
         this.path = file_path;
      }

      public String[] FileContents() throws IOException {
         FileReader fr = new FileReader(this.path);
         BufferedReader textReader = new BufferedReader(fr);
         Throwable var4 = null;

         String[] textData;
         try {
            int x2 = this.readLines();
            textData = new String[x2];

            for(int i = 0; i < x2; ++i) {
               textData[i] = textReader.readLine();
            }
         } catch (Throwable var14) {
            var4 = var14;
            throw var14;
         } finally {
            if(textReader != null) {
               if(var4 != null) {
                  try {
                     textReader.close();
                  } catch (Throwable var13) {
                     var4.addSuppressed(var13);
                  }
               } else {
                  textReader.close();
               }
            }

         }

         return textData;
      }

      public int readLines() throws IOException {
         FileReader file_to_read = new FileReader(this.path);
         BufferedReader bf = new BufferedReader(file_to_read);
         Throwable var4 = null;

         int numberOfLines;
         try {
            for(numberOfLines = 0; bf.readLine() != null; ++numberOfLines) {
               ;
            }
         } catch (Throwable var13) {
            var4 = var13;
            throw var13;
         } finally {
            if(bf != null) {
               if(var4 != null) {
                  try {
                     bf.close();
                  } catch (Throwable var12) {
                     var4.addSuppressed(var12);
                  }
               } else {
                  bf.close();
               }
            }

         }

         return numberOfLines;
      }
   }

   private static class FreqTable {
      private String path;
      private ArrayList loci_list;
      private LinkedHashMap freq_table;

      FreqTable(String file_path) {
         this.path = file_path;
         this.loci_list = new ArrayList();
         this.freq_table = new LinkedHashMap();

         try {
            NocitAlgorithm.OpenAndReadFile e = new NocitAlgorithm.OpenAndReadFile(this.path);
            String[] fileLines = e.FileContents();

            int i;
            String line;
            String[] parts;
            String locus;
            for(i = 1; i < fileLines.length; ++i) {
               line = fileLines[i];
               parts = line.split(",");
               locus = parts[0];
               if(!this.loci_list.contains(locus)) {
                  this.loci_list.add(locus);
               }
            }

            Iterator var14 = this.loci_list.iterator();

            while(var14.hasNext()) {
               line = (String)var14.next();
               this.freq_table.put(line, new LinkedHashMap());
            }

            for(i = 1; i < fileLines.length; ++i) {
               line = fileLines[i];
               parts = line.split(",");
               locus = parts[0];
               double a = Double.parseDouble(parts[1]);
               int allele = (int)Math.round(a * 10.0D);
               double frequency = Double.parseDouble(parts[2]);
               ((LinkedHashMap)this.freq_table.get(locus)).put(Integer.valueOf(allele), Double.valueOf(frequency));
            }
         } catch (IOException var13) {
            System.out.println(var13.getMessage());
         }

      }

      public ArrayList getFreqLociList() {
         return this.loci_list;
      }

      public LinkedHashMap getFreqTable() {
         return this.freq_table;
      }
   }

   private static class PosTwoExpFunction implements ParametricRealFunction {
      private PosTwoExpFunction() {
      }

      public double[] gradient(double x, double[] parameters) {
         double a = parameters[0];
         double b = parameters[1];
         double a_grad = Math.exp(b * x);
         double b_grad = a * x * Math.exp(b * x);
         double c_grad = 1.0D;
         if(parameters[2] < 0.0D) {
            c_grad = 0.0D;
         }

         double[] grad = new double[]{a_grad, b_grad, c_grad};
         return grad;
      }

      public double value(double x, double[] parameters) {
         return parameters[0] * Math.exp(parameters[1] * x) + Math.max(0.0D, parameters[2]);
      }

      // $FF: synthetic method
      PosTwoExpFunction(NocitAlgorithm.SyntheticClass_1 x0) {
         this();
      }
   }

   private static class TwoExpFunction implements ParametricRealFunction {
      private TwoExpFunction() {
      }

      public double[] gradient(double x, double[] parameters) {
         double a = parameters[0];
         double b = parameters[1];
         double a_grad = Math.exp(b * x);
         double b_grad = a * x * Math.exp(b * x);
         double c_grad = 1.0D;
         double[] grad = new double[]{a_grad, b_grad, c_grad};
         return grad;
      }

      public double value(double x, double[] parameters) {
         return parameters[0] * Math.exp(parameters[1] * x) + parameters[2];
      }

      // $FF: synthetic method
      TwoExpFunction(NocitAlgorithm.SyntheticClass_1 x0) {
         this();
      }
   }

   private static class ExpFunction implements ParametricRealFunction {
      private ExpFunction() {
      }

      public double[] gradient(double x, double[] parameters) {
         double a = parameters[0];
         double b = parameters[1];
         double a_grad = Math.exp(b * x);
         double b_grad = a * x * Math.exp(b * x);
         double[] grad = new double[]{a_grad, b_grad};
         return grad;
      }

      public double value(double x, double[] parameters) {
         double a = parameters[0];
         double b = parameters[1];
         return a * Math.exp(b * x);
      }

      // $FF: synthetic method
      ExpFunction(NocitAlgorithm.SyntheticClass_1 x0) {
         this();
      }
   }

   private static class LinFunction implements ParametricRealFunction {
      private LinFunction() {
      }

      public double[] gradient(double x, double[] parameters) {
         double b_grad = 1.0D;
         double[] grad = new double[]{x, b_grad};
         return grad;
      }

      public double value(double x, double[] parameters) {
         double a = parameters[0];
         double b = parameters[1];
         return a * x + b;
      }

      // $FF: synthetic method
      LinFunction(NocitAlgorithm.SyntheticClass_1 x0) {
         this();
      }
   }

   private static class CalibDataCollection {
      private String calib_path;
      private HashSet calibLoci;
      private HashSet masses;
      private LinkedHashMap freq_table;
      private HashMap trueCalibData;
      private HashMap noiseCalibData;
      private HashMap rstutterCalibData;
      private HashMap fstutterCalibData;
      private HashMap rstutterDOCalibData;
      private HashMap fstutterDOCalibData;
      private HashMap dropoutCalibData;

      CalibDataCollection(String calibpath, LinkedHashMap freqtable) {
         this.calib_path = calibpath;
         this.freq_table = freqtable;
         this.dropoutCalibData = new HashMap();
         this.rstutterDOCalibData = new HashMap();
         this.fstutterDOCalibData = new HashMap();
         this.rstutterCalibData = new HashMap();
         this.fstutterCalibData = new HashMap();
         this.noiseCalibData = new HashMap();
         this.trueCalibData = new HashMap();
         this.calibLoci = new HashSet();
         this.masses = new HashSet();
         HashMap trueValues = new HashMap();
         HashMap noiseValues = new HashMap();
         HashMap rstutterValues = new HashMap();
         HashMap fstutterValues = new HashMap();
         HashMap rstutterDO = new HashMap();
         HashMap fstutterDO = new HashMap();
         HashMap dropout = new HashMap();

         int val1;
         int rstutdo_count;
         int fstutdo_count;
         double val2;
         ArrayList var61;
         double var65;
         double var75;
         try {
            NocitAlgorithm.OpenAndReadFile i$ = new NocitAlgorithm.OpenAndReadFile(this.calib_path);
            String[] mass = i$.FileContents();

            int i;
            String i$1;
            String[] locus;
            for(i = 1; i < mass.length; ++i) {
               i$1 = mass[i];
               locus = i$1.split(",");
               String dovalues = locus[4];
               this.calibLoci.add(dovalues);
               double do_count = Double.parseDouble(locus[1]);
               this.masses.add(Double.valueOf(do_count));
            }

            Iterator var56 = this.masses.iterator();

            String var62;
            while(var56.hasNext()) {
               double var57 = ((Double)var56.next()).doubleValue();
               trueValues.put(Double.valueOf(var57), new HashMap());
               noiseValues.put(Double.valueOf(var57), new HashMap());
               rstutterValues.put(Double.valueOf(var57), new HashMap());
               fstutterValues.put(Double.valueOf(var57), new HashMap());
               dropout.put(Double.valueOf(var57), new HashMap());
               rstutterDO.put(Double.valueOf(var57), new HashMap());
               fstutterDO.put(Double.valueOf(var57), new HashMap());
               this.dropoutCalibData.put(Double.valueOf(var57), new HashMap());
               this.rstutterDOCalibData.put(Double.valueOf(var57), new HashMap());
               this.fstutterDOCalibData.put(Double.valueOf(var57), new HashMap());
               this.trueCalibData.put(Double.valueOf(var57), new HashMap());
               this.noiseCalibData.put(Double.valueOf(var57), new HashMap());
               this.rstutterCalibData.put(Double.valueOf(var57), new HashMap());
               this.fstutterCalibData.put(Double.valueOf(var57), new HashMap());
               Iterator var60 = this.calibLoci.iterator();

               while(var60.hasNext()) {
                  var62 = (String)var60.next();
                  ((HashMap)trueValues.get(Double.valueOf(var57))).put(var62, new ArrayList());
                  ((HashMap)noiseValues.get(Double.valueOf(var57))).put(var62, new ArrayList());
                  ((HashMap)rstutterValues.get(Double.valueOf(var57))).put(var62, new ArrayList());
                  ((HashMap)fstutterValues.get(Double.valueOf(var57))).put(var62, new ArrayList());
                  ((HashMap)dropout.get(Double.valueOf(var57))).put(var62, new ArrayList());
                  ((HashMap)rstutterDO.get(Double.valueOf(var57))).put(var62, new ArrayList());
                  ((HashMap)fstutterDO.get(Double.valueOf(var57))).put(var62, new ArrayList());
                  ((HashMap)this.trueCalibData.get(Double.valueOf(var57))).put(var62, new ArrayList());
                  ((HashMap)this.noiseCalibData.get(Double.valueOf(var57))).put(var62, new ArrayList());
                  ((HashMap)this.rstutterCalibData.get(Double.valueOf(var57))).put(var62, new ArrayList());
                  ((HashMap)this.fstutterCalibData.get(Double.valueOf(var57))).put(var62, new ArrayList());
               }
            }

            label397:
            for(i = 1; i < mass.length; ++i) {
               i$1 = mass[i];
               locus = i$1.split(",");
               var61 = new ArrayList();
               String[] var63 = locus;
               int do_obs = locus.length;

               for(int i$2 = 0; i$2 < do_obs; ++i$2) {
                  String final_do = var63[i$2];
                  if(final_do != null) {
                     var61.add(final_do);
                  }
               }

               var62 = (String)var61.get(4);
               var65 = Double.parseDouble((String)var61.get(1));
               ArrayList amel_alleles;
               ArrayList trueStats;
               ArrayList trueMean;
               ArrayList val;
               ArrayList trueStddev;
               HashSet trueAlleles;
               int index1;
               int fstutdovalues;
               ArrayList var66;
               int var87;
               if(!var62.contentEquals("AMEL")) {
                  var66 = new ArrayList();
                  amel_alleles = new ArrayList();
                  trueStats = new ArrayList();
                  trueMean = new ArrayList();
                  val = new ArrayList();
                  trueStddev = new ArrayList();
                  trueAlleles = new HashSet();
                  HashSet var73 = new HashSet();
                  var75 = Double.parseDouble((String)var61.get(2));
                  int var79 = (int)Math.round(var75 * 10.0D);
                  double var83 = Double.parseDouble((String)var61.get(3));
                  rstutdo_count = (int)Math.round(var83 * 10.0D);
                  if(var79 != rstutdo_count) {
                     trueAlleles.add(Integer.valueOf(var79));
                     trueAlleles.add(Integer.valueOf(rstutdo_count));

                     for(var87 = 6; var87 < var61.size(); var87 += 3) {
                        if(!((String)var61.get(var87)).equals("OL")) {
                           double var89 = Double.parseDouble((String)var61.get(var87));
                           index1 = (int)Math.round(var89 * 10.0D);
                           fstutdovalues = Integer.parseInt((String)var61.get(var87 + 2));
                           NocitAlgorithm.Peak var97 = new NocitAlgorithm.Peak(index1, fstutdovalues);
                           var66.add(var97);
                           amel_alleles.add(Integer.valueOf(var97.getAllele()));
                        }
                     }

                     Iterator var90 = var66.iterator();

                     while(true) {
                        int parent_height;
                        int rstutterStats;
                        int rstutterMean;
                        NocitAlgorithm.Peak var91;
                        Iterator var99;
                        NocitAlgorithm.Peak var100;
                        do {
                           if(!var90.hasNext()) {
                              var90 = var66.iterator();

                              while(var90.hasNext()) {
                                 var91 = (NocitAlgorithm.Peak)var90.next();
                                 if(!trueStddev.contains(var91)) {
                                    trueStats.add(var91);
                                    trueMean.add(Integer.valueOf(var91.getAllele()));
                                 }
                              }

                              var90 = trueAlleles.iterator();

                              label328:
                              while(true) {
                                 int var93;
                                 Iterator var102;
                                 NocitAlgorithm.Peak var104;
                                 do {
                                    label308:
                                    do {
                                       while(true) {
                                          while(var90.hasNext()) {
                                             var93 = ((Integer)var90.next()).intValue();
                                             var73.add(Integer.valueOf(var93));
                                             index1 = var93 - 10;
                                             fstutdovalues = var93 + 10;
                                             if(trueMean.contains(Integer.valueOf(var93))) {
                                                var99 = trueStats.iterator();

                                                while(var99.hasNext()) {
                                                   var100 = (NocitAlgorithm.Peak)var99.next();
                                                   if(var100.getAllele() == var93) {
                                                      parent_height = var100.getHeight();
                                                      ((ArrayList)((HashMap)trueValues.get(Double.valueOf(var65))).get(var62)).add(Integer.valueOf(parent_height));
                                                      ((ArrayList)((HashMap)dropout.get(Double.valueOf(var65))).get(var62)).add(Integer.valueOf(1));
                                                   }
                                                }

                                                if(!trueMean.contains(Integer.valueOf(index1))) {
                                                   ((ArrayList)((HashMap)rstutterDO.get(Double.valueOf(var65))).get(var62)).add(Integer.valueOf(0));
                                                } else {
                                                   var73.add(Integer.valueOf(index1));
                                                   if(!trueAlleles.contains(Integer.valueOf(index1)) && !trueAlleles.contains(Integer.valueOf(index1 - 10))) {
                                                      var99 = trueStats.iterator();

                                                      label298:
                                                      while(true) {
                                                         do {
                                                            if(!var99.hasNext()) {
                                                               break label298;
                                                            }

                                                            var100 = (NocitAlgorithm.Peak)var99.next();
                                                         } while(var100.getAllele() != var93);

                                                         parent_height = var100.getHeight();
                                                         var102 = trueStats.iterator();

                                                         while(var102.hasNext()) {
                                                            var104 = (NocitAlgorithm.Peak)var102.next();
                                                            rstutterStats = var104.getAllele();
                                                            if(rstutterStats == index1) {
                                                               rstutterMean = var104.getHeight();
                                                               val2 = (double)((float)rstutterMean / (float)parent_height);
                                                               ((ArrayList)((HashMap)rstutterValues.get(Double.valueOf(var65))).get(var62)).add(Double.valueOf(val2));
                                                               ((ArrayList)((HashMap)rstutterDO.get(Double.valueOf(var65))).get(var62)).add(Integer.valueOf(1));
                                                            }
                                                         }
                                                      }
                                                   }
                                                }

                                                if(trueMean.contains(Integer.valueOf(fstutdovalues))) {
                                                   var73.add(Integer.valueOf(fstutdovalues));
                                                   continue label308;
                                                }

                                                ((ArrayList)((HashMap)fstutterDO.get(Double.valueOf(var65))).get(var62)).add(Integer.valueOf(0));
                                             } else {
                                                ((ArrayList)((HashMap)dropout.get(Double.valueOf(var65))).get(var62)).add(Integer.valueOf(0));
                                             }
                                          }

                                          if(!this.freq_table.containsKey(var62)) {
                                             ((ArrayList)((HashMap)noiseValues.get(Double.valueOf(var65))).get(var62)).add(Integer.valueOf(0));
                                             continue label397;
                                          }

                                          Set var92 = ((LinkedHashMap)this.freq_table.get(var62)).keySet();
                                          Iterator var94 = trueStats.iterator();

                                          while(true) {
                                             do {
                                                do {
                                                   if(!var94.hasNext()) {
                                                      continue label397;
                                                   }

                                                   NocitAlgorithm.Peak var96 = (NocitAlgorithm.Peak)var94.next();
                                                   fstutdovalues = var96.getAllele();
                                                   fstutdo_count = var96.getHeight();
                                                } while(var73.contains(Integer.valueOf(fstutdovalues)));
                                             } while(!var92.contains(Integer.valueOf(fstutdovalues)) && !var92.contains(Integer.valueOf(fstutdovalues + 10)) && !var92.contains(Integer.valueOf(fstutdovalues - 10)));

                                             ((ArrayList)((HashMap)noiseValues.get(Double.valueOf(var65))).get(var62)).add(Integer.valueOf(fstutdo_count));
                                          }
                                       }
                                    } while(trueAlleles.contains(Integer.valueOf(fstutdovalues)));
                                 } while(trueAlleles.contains(Integer.valueOf(fstutdovalues + 10)));

                                 var99 = trueStats.iterator();

                                 while(true) {
                                    do {
                                       if(!var99.hasNext()) {
                                          continue label328;
                                       }

                                       var100 = (NocitAlgorithm.Peak)var99.next();
                                    } while(var100.getAllele() != var93);

                                    parent_height = var100.getHeight();
                                    var102 = trueStats.iterator();

                                    while(var102.hasNext()) {
                                       var104 = (NocitAlgorithm.Peak)var102.next();
                                       rstutterStats = var104.getAllele();
                                       if(rstutterStats == fstutdovalues) {
                                          rstutterMean = var104.getHeight();
                                          val2 = (double)((float)rstutterMean / (float)parent_height);
                                          ((ArrayList)((HashMap)fstutterValues.get(Double.valueOf(var65))).get(var62)).add(Double.valueOf(val2));
                                          ((ArrayList)((HashMap)fstutterDO.get(Double.valueOf(var65))).get(var62)).add(Integer.valueOf(1));
                                       }
                                    }
                                 }
                              }
                           }

                           var91 = (NocitAlgorithm.Peak)var90.next();
                           index1 = var91.getAllele();
                           fstutdovalues = Collections.frequency(amel_alleles, Integer.valueOf(index1));
                        } while(fstutdovalues <= 1);

                        var99 = var66.iterator();

                        while(var99.hasNext()) {
                           var100 = (NocitAlgorithm.Peak)var99.next();
                           parent_height = var100.getAllele();
                           if(index1 == parent_height) {
                              int fstut_do = var66.indexOf(var91);
                              int peakobj21 = var66.indexOf(var100);
                              if(fstut_do != peakobj21) {
                                 rstutterStats = var91.getHeight();
                                 rstutterMean = var100.getHeight();
                                 if(rstutterStats >= rstutterMean) {
                                    if(!val.contains(Integer.valueOf(parent_height))) {
                                       trueStddev.add(var100);
                                       val.add(Integer.valueOf(parent_height));
                                    }
                                 } else if(!val.contains(Integer.valueOf(index1))) {
                                    trueStddev.add(var91);
                                    val.add(Integer.valueOf(index1));
                                 }
                              }
                           }
                        }
                     }
                  }
               } else {
                  var66 = new ArrayList();
                  amel_alleles = new ArrayList();
                  trueStats = new ArrayList();
                  trueMean = new ArrayList();
                  val = new ArrayList();
                  trueStddev = new ArrayList();
                  trueAlleles = new HashSet();
                  String noiseStats = (String)var61.get(2);
                  String noiseMean = (String)var61.get(3);
                  trueAlleles.add(noiseStats);
                  trueAlleles.add(noiseMean);

                  for(val1 = 6; val1 < var61.size(); val1 += 3) {
                     if(!((String)var61.get(val1)).equals("OL")) {
                        String noiseStddev = (String)var61.get(val1);
                        int obj = Integer.parseInt((String)var61.get(val1 + 2));
                        NocitAlgorithm.AMEL_Peak rstutdovalues = new NocitAlgorithm.AMEL_Peak(noiseStddev, obj);
                        var66.add(rstutdovalues);
                        amel_alleles.add(rstutdovalues.getAllele());
                     }
                  }

                  Iterator var77 = var66.iterator();

                  while(true) {
                     NocitAlgorithm.AMEL_Peak var78;
                     String var82;
                     do {
                        if(!var77.hasNext()) {
                           var77 = var66.iterator();

                           while(var77.hasNext()) {
                              var78 = (NocitAlgorithm.AMEL_Peak)var77.next();
                              if(!trueStddev.contains(var78)) {
                                 trueStats.add(var78);
                                 trueMean.add(var78.getAllele());
                              }
                           }

                           if(noiseStats.contentEquals(noiseMean)) {
                              if(!trueMean.contains("Y") || trueAlleles.contains("Y")) {
                                 continue label397;
                              }

                              var77 = trueStats.iterator();

                              while(true) {
                                 if(!var77.hasNext()) {
                                    continue label397;
                                 }

                                 var78 = (NocitAlgorithm.AMEL_Peak)var77.next();
                                 if(var78.getAllele().contentEquals("Y")) {
                                    ((ArrayList)((HashMap)noiseValues.get(Double.valueOf(var65))).get(var62)).add(Integer.valueOf(var78.getHeight()));
                                 }
                              }
                           }

                           var77 = trueAlleles.iterator();

                           while(true) {
                              while(true) {
                                 if(!var77.hasNext()) {
                                    continue label397;
                                 }

                                 String var81 = (String)var77.next();
                                 if(trueMean.contains(var81)) {
                                    Iterator var84 = trueStats.iterator();

                                    while(var84.hasNext()) {
                                       NocitAlgorithm.AMEL_Peak var86 = (NocitAlgorithm.AMEL_Peak)var84.next();
                                       if(var81.contentEquals(var86.getAllele())) {
                                          var87 = var86.getHeight();
                                          ((ArrayList)((HashMap)trueValues.get(Double.valueOf(var65))).get(var62)).add(Integer.valueOf(var87));
                                          ((ArrayList)((HashMap)dropout.get(Double.valueOf(var65))).get(var62)).add(Integer.valueOf(1));
                                       }
                                    }
                                 } else {
                                    ((ArrayList)((HashMap)dropout.get(Double.valueOf(var65))).get(var62)).add(Integer.valueOf(0));
                                 }
                              }
                           }
                        }

                        var78 = (NocitAlgorithm.AMEL_Peak)var77.next();
                        var82 = var78.getAllele();
                        rstutdo_count = Collections.frequency(amel_alleles, var82);
                     } while(rstutdo_count <= 1);

                     Iterator rstut_obs = var66.iterator();

                     while(rstut_obs.hasNext()) {
                        NocitAlgorithm.AMEL_Peak peakobj2 = (NocitAlgorithm.AMEL_Peak)rstut_obs.next();
                        String rstut_do = peakobj2.getAllele();
                        if(var82.contentEquals(rstut_do)) {
                           index1 = var66.indexOf(var78);
                           fstutdovalues = var66.indexOf(peakobj2);
                           if(index1 != fstutdovalues) {
                              fstutdo_count = var78.getHeight();
                              int fstut_obs = peakobj2.getHeight();
                              if(fstutdo_count >= fstut_obs) {
                                 if(!val.contains(rstut_do)) {
                                    trueStddev.add(peakobj2);
                                    val.add(rstut_do);
                                 }
                              } else if(!val.contains(var82)) {
                                 trueStddev.add(var78);
                                 val.add(var82);
                              }
                           }
                        }
                     }
                  }
               }
            }
         } catch (IOException var53) {
            System.out.println(var53.getMessage());
         }

         Iterator var54 = this.masses.iterator();

         label225:
         while(var54.hasNext()) {
            double var55 = ((Double)var54.next()).doubleValue();
            Iterator var58 = this.calibLoci.iterator();

            while(true) {
               String var59;
               do {
                  if(!var58.hasNext()) {
                     continue label225;
                  }

                  var59 = (String)var58.next();
                  var61 = (ArrayList)((HashMap)dropout.get(Double.valueOf(var55))).get(var59);
                  int var64 = Collections.frequency(var61, Integer.valueOf(0));
                  var65 = (double)var61.size();
                  double var67 = (double)var64 / var65;
                  ((HashMap)this.dropoutCalibData.get(Double.valueOf(var55))).put(var59, Double.valueOf(var67));
                  DescriptiveStatistics var68 = new DescriptiveStatistics();
                  Iterator var69 = ((ArrayList)((HashMap)trueValues.get(Double.valueOf(var55))).get(var59)).iterator();

                  while(var69.hasNext()) {
                     int var71 = ((Integer)var69.next()).intValue();
                     var68.addValue((double)var71);
                  }

                  double var70 = var68.getMean();
                  double var72 = var68.getStandardDeviation();
                  ((ArrayList)((HashMap)this.trueCalibData.get(Double.valueOf(var55))).get(var59)).add(Double.valueOf(var70));
                  ((ArrayList)((HashMap)this.trueCalibData.get(Double.valueOf(var55))).get(var59)).add(Double.valueOf(var72));
                  DescriptiveStatistics var74 = new DescriptiveStatistics();
                  Iterator var76 = ((ArrayList)((HashMap)noiseValues.get(Double.valueOf(var55))).get(var59)).iterator();

                  while(var76.hasNext()) {
                     val1 = ((Integer)var76.next()).intValue();
                     var74.addValue((double)val1);
                  }

                  var75 = var74.getMean();
                  double var80 = var74.getStandardDeviation();
                  ((ArrayList)((HashMap)this.noiseCalibData.get(Double.valueOf(var55))).get(var59)).add(Double.valueOf(var75));
                  ((ArrayList)((HashMap)this.noiseCalibData.get(Double.valueOf(var55))).get(var59)).add(Double.valueOf(var80));
               } while(var59.contentEquals("AMEL"));

               ArrayList var85 = (ArrayList)((HashMap)rstutterDO.get(Double.valueOf(var55))).get(var59);
               rstutdo_count = Collections.frequency(var85, Integer.valueOf(0));
               double var88 = (double)var85.size();
               double var95 = (double)rstutdo_count / var88;
               ((HashMap)this.rstutterDOCalibData.get(Double.valueOf(var55))).put(var59, Double.valueOf(var95));
               ArrayList var98 = (ArrayList)((HashMap)fstutterDO.get(Double.valueOf(var55))).get(var59);
               fstutdo_count = Collections.frequency(var98, Integer.valueOf(0));
               double var101 = (double)var98.size();
               double var103 = (double)fstutdo_count / var101;
               ((HashMap)this.fstutterDOCalibData.get(Double.valueOf(var55))).put(var59, Double.valueOf(var103));
               DescriptiveStatistics var105 = new DescriptiveStatistics();
               Iterator var106 = ((ArrayList)((HashMap)rstutterValues.get(Double.valueOf(var55))).get(var59)).iterator();

               while(var106.hasNext()) {
                  val2 = ((Double)var106.next()).doubleValue();
                  var105.addValue(val2);
               }

               double var107 = var105.getMean();
               double rstutterStddev = var105.getStandardDeviation();
               ((ArrayList)((HashMap)this.rstutterCalibData.get(Double.valueOf(var55))).get(var59)).add(Double.valueOf(var107));
               ((ArrayList)((HashMap)this.rstutterCalibData.get(Double.valueOf(var55))).get(var59)).add(Double.valueOf(rstutterStddev));
               DescriptiveStatistics fstutterStats = new DescriptiveStatistics();
               Iterator fstutterMean = ((ArrayList)((HashMap)fstutterValues.get(Double.valueOf(var55))).get(var59)).iterator();

               while(fstutterMean.hasNext()) {
                  double val3 = ((Double)fstutterMean.next()).doubleValue();
                  fstutterStats.addValue(val3);
               }

               double var108 = fstutterStats.getMean();
               double fstutterStddev = fstutterStats.getStandardDeviation();
               ((ArrayList)((HashMap)this.fstutterCalibData.get(Double.valueOf(var55))).get(var59)).add(Double.valueOf(var108));
               ((ArrayList)((HashMap)this.fstutterCalibData.get(Double.valueOf(var55))).get(var59)).add(Double.valueOf(fstutterStddev));
            }
         }

      }

      public HashMap getDOCalibData() {
         return this.dropoutCalibData;
      }

      public HashMap getfstutDOCalibData() {
         return this.fstutterDOCalibData;
      }

      public HashMap getfstutterCalibData() {
         return this.fstutterCalibData;
      }

      public HashSet getLoci() {
         return this.calibLoci;
      }

      public HashSet getMasses() {
         return this.masses;
      }

      public HashMap getNoiseCalibData() {
         return this.noiseCalibData;
      }

      public HashMap getrstutDOCalibData() {
         return this.rstutterDOCalibData;
      }

      public HashMap getrstutterCalibData() {
         return this.rstutterCalibData;
      }

      public HashMap getTrueCalibData() {
         return this.trueCalibData;
      }
   }

   private static class Callable_amel_object {
      private double thread_value;
      private double max_val;
      private Set max_alls;

      Callable_amel_object(double threadval, double maxval, Set maxalls) {
         this.thread_value = threadval;
         this.max_val = maxval;
         this.max_alls = maxalls;
      }

      public Set getmaxalls() {
         return this.max_alls;
      }

      public double getmaxvalue() {
         return this.max_val;
      }

      public double getvalue() {
         return this.thread_value;
      }
   }

   private static class Callable_loci_object {
      private double thread_value;
      private double max_val;
      private Set max_alls;

      Callable_loci_object(double threadval, double maxval, Set maxalls) {
         this.thread_value = threadval;
         this.max_val = maxval;
         this.max_alls = maxalls;
      }

      public Set getmaxalls() {
         return this.max_alls;
      }

      public double getmaxvalue() {
         return this.max_val;
      }

      public double getvalue() {
         return this.thread_value;
      }
   }

   private static class AMEL_Peak {
      private String allele;
      private int height;

      public AMEL_Peak(String assignallele, int assignheight) {
         this.allele = assignallele;
         this.height = assignheight;
      }

      public String getAllele() {
         return this.allele;
      }

      public int getHeight() {
         return this.height;
      }

      public String toString() {
         return this.allele + "," + this.height;
      }
   }

   private static class Peak {
      private int allele;
      private int height;

      Peak(int assignallele, int assignheight) {
         this.allele = assignallele;
         this.height = assignheight;
      }

      public int getAllele() {
         return this.allele;
      }

      public int getHeight() {
         return this.height;
      }

      public String toString() {
         return this.allele + "," + this.height;
      }
   }
}
