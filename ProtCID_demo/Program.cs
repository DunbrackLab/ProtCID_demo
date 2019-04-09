using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Threading.Tasks;
using ProtCidSettingsLib;

namespace ProtCID_demo
{
    class Program
    {
        /// <summary>
        /// input
        /// 1. -infile a text file containing a list of PDB, one entry per line
        /// 2. -datadir a directory name where all files can be saved
        /// 3. -alnfile a user alignment file to provide the correspondence of residue numbers
        /// 4. -groupname a name for coordinate files of clusters
        /// </summary>
        /// <param name="args"></param>
        static void Main(string[] args)
        {
            string exeDir = System.Reflection.Assembly.GetExecutingAssembly().GetName().CodeBase;
            ProtCidSettings.applicationStartPath = Path.GetDirectoryName(exeDir).Replace("file:\\", "");
            ProtCidSettings.paramFile = Path.Combine(ProtCidSettings.applicationStartPath, "Settings\\parameters.xml");
            ProtCidSettings.dirFile = Path.Combine(ProtCidSettings.applicationStartPath, "Settings\\dirSettings.xml");
            ProtCidSettings.symOpsFile = Path.Combine(ProtCidSettings.applicationStartPath, "Settings\\symOps.xml");
            ProtCidSettings.crystMethodFile = Path.Combine(ProtCidSettings.applicationStartPath, "Settings\\CrystMethods.txt");

            // AppDomain.CurrentDomain.BaseDirectory;            
            // default directory settings and file names 
            string dataDir = @"demo_data";
            string entryFile = @"demo_data\ls-pdb_ST1A1.txt";
            //          string entryFile = @"demo_data\ls-pdb_RAS.txt";
            string alignFile = @"demo_data\RasMonomers_clustalO.aln";
            string outGroupName = "sulf";

            bool hasAlignFile = false;

            for (int i = 0; i < args.Length; i += 2)
            {
                switch (args[i].ToLower())
                {
                    case "-infile":
                        entryFile = args[i + 1];
                        break;

                    case "-datadir":
                        dataDir = args[i + 1];
                        break;

                    case "-alnfile":
                        alignFile = args[i + 1];
                        hasAlignFile = true;
                        break;

                    case "-groupname":
                        outGroupName = args[i + 1];
                        break;

                    default:
                        break;
                }
            }

            ProtCidSettings.tempDir = Path.Combine(dataDir, "xtal_temp");
            if (!Directory.Exists(ProtCidSettings.tempDir))
            {
                Directory.CreateDirectory(ProtCidSettings.tempDir);
            }

            string logFile = Path.Combine(dataDir, "log.txt");
            ProtCidSettings.logWriter = new StreamWriter(logFile, true);

            if (!File.Exists(entryFile))
            {
                Console.WriteLine("The entry file : " + entryFile + " is not exist. " +
                    " You must provide a valid text file containing a list of PDBs, one PDB per line");
            }
            else
            {
                //            string[] pdbIds = { "1ls6", "1z28", "2a3r", "3u3r", "3u3o", "3qvu", "4gra" };
                InterfaceClustering crystInterfaceCluster = new InterfaceClustering(dataDir, outGroupName);
                if (hasAlignFile)
                {
                    crystInterfaceCluster.DemonstrateProtCidMainFunctions(entryFile, alignFile);
                }
                else
                {
                    crystInterfaceCluster.DemonstrateProtCidMainFunctions(entryFile);
                }
            }
            try
            {
                Directory.Delete(ProtCidSettings.tempDir, true);
            }
            catch { }
            ProtCidSettings.logWriter.Close();

            Console.WriteLine("Done!");
            Console.ReadLine();
        }
    }
}
