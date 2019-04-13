using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Data;
using System.IO;
using CrystalInterfaceLib.Contacts;
using CrystalInterfaceLib.Crystal;
using CrystalInterfaceLib.StructureComp;
using InterfaceClusterLib.Clustering;
using InterfaceClusterLib.PymolScript;
using InterfaceClusterLib.AuxFuncs;
using AuxFuncLib;
using BuQueryLib;
using ProtCidSettingsLib;

namespace ProtCID_demo
{
    public class InterfaceClustering
    {
        #region variables
        private string dataDir = "";
        private string outGroupName = "";
        private InterfacesComp interfaceComp = new InterfacesComp();
        private InterfaceCluster interfaceCluster = new InterfaceCluster();
        private BiolUnitQuery buQuery = new BiolUnitQuery();
        private double entryHomoQCutoff = 0.9;
        private InterfaceAlignPymolScript pmlScriptWriter = new InterfaceAlignPymolScript();
        private FileCompress fileCompress = new FileCompress();
        private string residueNumberingType = "auth";

        private string interfaceFileDir = "";
        private string xmlFileDir = "";
        private string textFileDir = "";
        private string clusterCoordFileDir = "";

        public InterfaceClustering (string inDataDir, string outName)
        {
            dataDir = inDataDir;
            outGroupName = outName;
            InitializeDirectorySetting();
        }

        public InterfaceClustering()
        {
            dataDir = "";   // the directory in executable program
            InitializeDirectorySetting();
        }

        /// <summary>
        /// 
        /// </summary>
        private void InitializeDirectorySetting ()
        {
            if (outGroupName != "")
            {
                dataDir = Path.Combine(dataDir, outGroupName);               
            }
            interfaceFileDir = Path.Combine(dataDir, "CrystInterfaces");
            if (!Directory.Exists(interfaceFileDir))
            {
                Directory.CreateDirectory(interfaceFileDir);
            }
            xmlFileDir = Path.Combine(dataDir, "xml");
            if (!Directory.Exists(xmlFileDir))
            {
                Directory.CreateDirectory(xmlFileDir);
            }
            textFileDir = Path.Combine(dataDir, "textfiles");
            if (!Directory.Exists(textFileDir))
            {
                Directory.CreateDirectory(textFileDir);
            }
            clusterCoordFileDir = Path.Combine(dataDir, "clusterCoordinates");
            if (!Directory.Exists(clusterCoordFileDir))
            {
                Directory.CreateDirectory(clusterCoordFileDir);
            }
        }
        #endregion

        #region main interfaces of protcid demo
        /// <summary>
        /// 
        /// </summary>
        /// <param name="entryLsFile">text file containing a list of pdb, one entry per line</param>
        /// <param name="dataDir">the directory where all files to be saved</param>
        public void DemonstrateProtCidMainFunctions (string entryLsFile)
        {
            string[] pdbIds = ReadEntries(entryLsFile);

            DemonstrateProtCidMainFunctions(pdbIds);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="entryLsFile">text file containing a list of pdb, one entry per line</param>
        /// <param name="seqAlignmentFile"></param>
        /// <param name="dataDir">the directory where all files to be saved</param>
        public void DemonstrateProtCidMainFunctions(string entryLsFile, string seqAlignmentFile)
        {
            string[] pdbIds = ReadEntries(entryLsFile);
            DataTable seqAlignmentTable = ReadAlignmentFile(seqAlignmentFile, pdbIds);

            DemonstrateProtCidMainFunctions(pdbIds, seqAlignmentTable);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="pdbIds"></param>
        /// <param name="dataDir"></param>
        public void DemonstrateProtCidMainFunctions (string[] pdbIds)
        {
            residueNumberingType = "auth";

            Console.WriteLine("Calculate interfaces from crystals");
            DataTable[] asuTables = null;
            // compute interfaces in PDB crystals and generate interface coordinate files
            Dictionary<string, InterfaceChains[]> entryCrystInterfacesDict = CalculateCrystInterfaces(pdbIds, xmlFileDir, interfaceFileDir, out asuTables);

            Console.WriteLine("Output details of input structures");
            PrintEntryAsuInfoToTextFiles(asuTables, textFileDir);

            Console.WriteLine("Output the details of interfaces from crystals");
            string interfaceTextFile = Path.Combine(textFileDir, "CrystInterfaces.txt");
            // print crystal interfaces to a text file
            PrintEntryCrystInterfacesToFile(entryCrystInterfacesDict, interfaceTextFile, asuTables[0]);

            Console.WriteLine("Retreive unique interfaces");
            DataTable entryInterfaceCompTable = null;
            GetEntryUniqueInterfaces(entryCrystInterfacesDict, out entryInterfaceCompTable);

            Console.WriteLine("Output similarity Q scores of same-entry interfaces to a text file");
            string sameEntryInterfaceCompFile = Path.Combine(textFileDir, "SameEntryInterfaceCompInfo.txt");
            PrintSameEntryInterfaceCompToTextFile(entryInterfaceCompTable, sameEntryInterfaceCompFile);

            Console.WriteLine("Representative interfaces of each entry and its redundant interfaces");
            // similar interfaces for representative interfaces 
            Dictionary<string, string[]> entryRepHomoInterfacesDict = GetEntryHomoInterfacesDict(entryCrystInterfacesDict, entryInterfaceCompTable);

            Console.WriteLine("Calculate similarity Q scores of interfaces between entries");
            // Calculate similarity Q scores between two entries
            DataTable difEntryInterfaceCompTable = CompareDiffEntryInterfaces(entryCrystInterfacesDict);

            Console.WriteLine("Output similarity Q scores of diff-entry interfaces to a text file");
            string diffEntryInterfaceCompFile = Path.Combine(textFileDir, "DiffEntryInterfaceCompInfo.txt");
            PrintDiffEntryInterfaceCompToTextFile(difEntryInterfaceCompTable, diffEntryInterfaceCompFile);

            Console.WriteLine("entry and spacegroup-asu");
            // group structures by crystal form (space group + ASU)
            Dictionary<string, string> entryCfGroupDict = GetEntrySpacegroupAsuDict(pdbIds, asuTables[0], asuTables[2]);

            Console.WriteLine("Clustering interfaces");
            // clustering interfaces based on similarity Q scores
            DataTable clusterTable = ClusterInterfaces(difEntryInterfaceCompTable, entryInterfaceCompTable, entryCfGroupDict, entryRepHomoInterfacesDict);
            Console.WriteLine("Output interface clusters to a text file");
            // output clusters to a text file
            string clusterTextFile = Path.Combine(textFileDir, "InterfaceClusters.txt");
            PrintClusterTableToTextFile(clusterTable, clusterTextFile);

            Console.WriteLine("Compile coordinates and pymol scripts of clusters");
            CompileClusterCoordinates(outGroupName, clusterTable, clusterCoordFileDir, interfaceFileDir);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="pdbIds"></param>
        /// <param name="dataDir"></param>
        private void DemonstrateProtCidMainFunctions(string[] pdbIds, DataTable seqAlignmentTable)
        {
            residueNumberingType = "xml"; // the sequential residue numbering

            Console.WriteLine("Calculate interfaces from crystals");
            DataTable[] asuTables = null;
            // compute interfaces in PDB crystals and generate interface coordinate files
            Dictionary<string, InterfaceChains[]> entryCrystInterfacesDict = CalculateCrystInterfaces(pdbIds, xmlFileDir, interfaceFileDir, out asuTables);

            Console.WriteLine("Output details of input structures");
            PrintEntryAsuInfoToTextFiles(asuTables, textFileDir);

            FillAlignmentTableEntityChainIds (seqAlignmentTable, asuTables[0]);

            Console.WriteLine("Output the details of interfaces from crystals");
            string interfaceTextFile = Path.Combine(textFileDir, "CrystInterfaces.txt");
            // print crystal interfaces to a text file
            PrintEntryCrystInterfacesToFile(entryCrystInterfacesDict, interfaceTextFile, asuTables[0]);

            Console.WriteLine("Retreive unique interfaces");
            DataTable entryInterfaceCompTable = null;
            GetEntryUniqueInterfaces(entryCrystInterfacesDict, out entryInterfaceCompTable);

            Console.WriteLine("Output similarity Q scores of same-entry interfaces to a text file");
            string sameEntryInterfaceCompFile = Path.Combine(textFileDir, "SameEntryInterfaceCompInfo.txt");
            PrintSameEntryInterfaceCompToTextFile(entryInterfaceCompTable, sameEntryInterfaceCompFile);

            Console.WriteLine("Representative interfaces of each entry and its redundant interfaces");
            // similar interfaces for representative interfaces 
            Dictionary<string, string[]> entryRepHomoInterfacesDict = GetEntryHomoInterfacesDict(entryCrystInterfacesDict, entryInterfaceCompTable);

            Console.WriteLine("Calculate similarity Q scores of interfaces between entries");
            // Calculate similarity Q scores between two entries
            DataTable difEntryInterfaceCompTable = CompareDiffEntryInterfaces(entryCrystInterfacesDict, seqAlignmentTable);

            Console.WriteLine("Output similarity Q scores of diff-entry interfaces to a text file");
            string diffEntryInterfaceCompFile = Path.Combine(textFileDir, "DiffEntryInterfaceCompInfo.txt");
            PrintDiffEntryInterfaceCompToTextFile(difEntryInterfaceCompTable, diffEntryInterfaceCompFile);

            Console.WriteLine("entry and spacegroup-asu");
            // group structures by crystal form (space group + ASU)
            Dictionary<string, string> entryCfGroupDict = GetEntrySpacegroupAsuDict(pdbIds, asuTables[0], asuTables[2]);

            Console.WriteLine("Clustering interfaces");
            // clustering interfaces based on similarity Q scores
            DataTable clusterTable = ClusterInterfaces(difEntryInterfaceCompTable, entryInterfaceCompTable, entryCfGroupDict, entryRepHomoInterfacesDict);
            Console.WriteLine("Output interface clusters to a text file");
            // output clusters to a text file
            string clusterTextFile = Path.Combine(textFileDir, "InterfaceClusters.txt");
            PrintClusterTableToTextFile(clusterTable, clusterTextFile);

            Console.WriteLine("Compile coordinates and pymol scripts of clusters");
            CompileClusterCoordinates(outGroupName, clusterTable, clusterCoordFileDir, interfaceFileDir);
        }
        #endregion

        #region interfaces in crystals
        /// <summary>
        /// 
        /// </summary>
        /// <param name="pdbIds"></param>
        /// <param name="dataDir"></param>
        /// <param name="interfaceFileDir"></param>
        /// <returns></returns>
        public Dictionary<string, InterfaceChains[]> CalculateCrystInterfaces(string[] pdbIds, string xmlFileDir, string interfaceFileDir, out DataTable[] asuTables)
        {
            CrystInterfacesGen interfaceGen = new CrystInterfacesGen();
            XmlFileParser xmlParser = new XmlFileParser(xmlFileDir);
            Dictionary<string, InterfaceChains[]> entryCrystInterfacesDict = new Dictionary<string, InterfaceChains[]>();
            string pdbXmlFile = "";
            asuTables = new DataTable[3];
            foreach (string pdbId in pdbIds)
            {
                Console.WriteLine(pdbId);
                pdbXmlFile = Path.Combine(xmlFileDir, pdbId + ".xml.gz");
                try
                {
                    DataTable[] entryAsuTables = xmlParser.RetrieveEntryInfoFromXmlFile(pdbXmlFile);
                    ParseHelper.AddNewTableToExistTable(entryAsuTables[0], ref asuTables[0]); // asymmetric unit: entityid, asymid, authorchain
                    ParseHelper.AddNewTableToExistTable(entryAsuTables[1], ref asuTables[1]);  // db reference
                    ParseHelper.AddNewTableToExistTable(entryAsuTables[2], ref asuTables[2]);  // entry info: method, rfactor, space group
                    InterfaceChains[] crystInterfaces = interfaceGen.GenerateInterfacesFromCryst(pdbId, xmlFileDir, interfaceFileDir, entryAsuTables[0], residueNumberingType);
                    entryCrystInterfacesDict.Add(pdbId, crystInterfaces);
                }
                catch (Exception ex)
                {
                    Console.WriteLine("Calculating interfaces in the crystal " + pdbId + " error: " + ex.Message);
                }
            }
            return entryCrystInterfacesDict;
        }
       
        /// <summary>
        /// 
        /// </summary>
        /// <param name="entryCrystInterfacesDict"></param>
        /// <param name="entryInterfaceCompTable"></param>
        public void GetEntryUniqueInterfaces (Dictionary<string, InterfaceChains[]> entryCrystInterfacesDict, out DataTable entryInterfaceCompTable)
        {
      //      Dictionary<string, InterfaceChains[]> entryUniqueInterfaceDict = new Dictionary<string, InterfaceChains[]>();
            entryInterfaceCompTable = InitializeSameEntryInterfaceCompTable();
            foreach (string pdbId in entryCrystInterfacesDict.Keys)
            {
                InterfaceChains[] entryInterfaces = entryCrystInterfacesDict[pdbId];
                InterfacePairInfo[] compInfo = interfaceComp.CompareInterfacesWithinCrystal(ref entryInterfaces);
                foreach (InterfacePairInfo interfacePair in compInfo)
                {
                    DataRow dataRow = entryInterfaceCompTable.NewRow();
                    dataRow["PdbID"] = pdbId;
                    dataRow["InterfaceID1"] = interfacePair.interfaceInfo1.interfaceId;
                    dataRow["InterfaceID2"] = interfacePair.interfaceInfo2.interfaceId;
                    dataRow["Qscore"] = interfacePair.qScore;
                    entryInterfaceCompTable.Rows.Add(dataRow);
                }
            }
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="entryUniqueInterfacesDict"></param>
        /// <param name="entryInterfaceCompTable"></param>
        /// <returns></returns>
        private Dictionary<string, string[]> GetEntryHomoInterfacesDict (Dictionary<string, InterfaceChains[]> entryUniqueInterfacesDict, DataTable entryInterfaceCompTable)
        {
            Dictionary<string, string[]> entryHomoInterfacesDict = new Dictionary<string, string[]>();
            List<int> entryRepInterfaceIdList = new List<int>();
            List<int> entryInterfaceIdList = new List<int>();
            foreach (string pdbId in entryUniqueInterfacesDict.Keys)
            {
                entryRepInterfaceIdList.Clear();
                entryInterfaceIdList.Clear();
                foreach (InterfaceChains crystInterface in entryUniqueInterfacesDict[pdbId])
                {
                    entryRepInterfaceIdList.Add(crystInterface.interfaceId);
                }
                entryInterfaceIdList.AddRange(entryRepInterfaceIdList);
                foreach (int repInterfaceId in entryRepInterfaceIdList)
                {
                    int[] homoInterfaceIds = GetHomoInterfaceIds(pdbId, repInterfaceId, entryInterfaceCompTable);
                    List<string> homoInterfaceList = new List<string>();
                    foreach (int homoInterfaceId in homoInterfaceIds)
                    {
                        if (!entryInterfaceIdList.Contains(homoInterfaceId))
                        {
                            homoInterfaceList.Add(pdbId + "_" + homoInterfaceId);
                            entryInterfaceIdList.Add(homoInterfaceId);
                        }
                    }
                    if (homoInterfaceList.Count > 0)
                    {
                        entryHomoInterfacesDict.Add(pdbId + "_" + repInterfaceId, homoInterfaceList.ToArray());
                    }
                }
            }
            return entryHomoInterfacesDict;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="pdbId"></param>
        /// <param name="repInterfaceId"></param>
        /// <param name="entryInterfaceCompTable"></param>
        /// <returns></returns>
        private int[] GetHomoInterfaceIds (string pdbId, int repInterfaceId,  DataTable entryInterfaceCompTable)
        {
            DataRow[] simInterfaceRows = entryInterfaceCompTable.Select(string.Format ("PdbID = '{0}' AND InterfaceID1 = {1} AND Qscore >= {2}", pdbId, repInterfaceId, entryHomoQCutoff));
            List<int> homoInterfaceIdList = new List<int>();
            int interfaceId = 0;
            foreach (DataRow compRow in simInterfaceRows)
            {
                interfaceId = Convert.ToInt32(compRow["InterfaceID2"].ToString ());
                if (! homoInterfaceIdList.Contains (interfaceId))
                {
                    homoInterfaceIdList.Add(interfaceId);
                }
            }

            simInterfaceRows = entryInterfaceCompTable.Select(string.Format ("PdbID = '{0}' AND InterfaceID2 = {1} AND Qscore >= {2}", pdbId, repInterfaceId, entryHomoQCutoff));
            foreach (DataRow compRow in simInterfaceRows)
            {
                interfaceId = Convert.ToInt32(compRow["InterfaceID1"].ToString ());
                if (! homoInterfaceIdList.Contains (interfaceId))
                {
                    homoInterfaceIdList.Add(interfaceId);
                }
            }

            return homoInterfaceIdList.ToArray();
        }
        #endregion

        #region initialize data tables
        /// <summary>
        /// Data table for Q scores between interfaces of one entry
        /// </summary>
        /// <returns></returns>
        private DataTable InitializeSameEntryInterfaceCompTable ()
        {
            DataTable sameEntryInterfaceCompTable = new DataTable("SameEntryInterfaceComp");
   //         string[] columns = {"PdbID", "InterfaceID1", "InterfaceID2", "Qscore"};
            sameEntryInterfaceCompTable.Columns.Add(new DataColumn("PdbID"));
            sameEntryInterfaceCompTable.Columns.Add(new DataColumn("InterfaceID1", Type.GetType ("System.Int32")));
            sameEntryInterfaceCompTable.Columns.Add(new DataColumn("InterfaceID2", Type.GetType("System.Int32")));
            sameEntryInterfaceCompTable.Columns.Add(new DataColumn("Qscore", Type.GetType("System.Double")));            
            return sameEntryInterfaceCompTable;
        }

        /// <summary>
        /// Data table for Q scores between interfaces of two entries
        /// </summary>
        /// <returns></returns>
        private DataTable InitializeDifEntryInterfaceCompTable()
        {
            DataTable difEntryInterfaceCompTable = new DataTable("DiffEntryInterfaceComp");
            string[] columns = { "PdbID1", "InterfaceID1", "PdbID2", "InterfaceID2", "Qscore"};
            difEntryInterfaceCompTable.Columns.Add(new DataColumn("PdbID1"));
            difEntryInterfaceCompTable.Columns.Add(new DataColumn("InterfaceID1", Type.GetType("System.Int32")));
            difEntryInterfaceCompTable.Columns.Add(new DataColumn("PdbID2"));
            difEntryInterfaceCompTable.Columns.Add(new DataColumn("InterfaceID2", Type.GetType("System.Int32")));
            difEntryInterfaceCompTable.Columns.Add(new DataColumn("Qscore", Type.GetType("System.Double")));
            return difEntryInterfaceCompTable;
        }
        #endregion

        #region similarity scores
        /// <summary>
        /// 
        /// </summary>
        /// <param name="entryCrystInterfacesDict"></param>
        /// <returns></returns>
        private DataTable CompareDiffEntryInterfaces(Dictionary<string, InterfaceChains[]> entryCrystInterfacesDict)
        {
            DataTable difEntryInterfaceCompTable = InitializeDifEntryInterfaceCompTable();
            List<string> pdbIdList = new List<string>(entryCrystInterfacesDict.Keys);
            pdbIdList.Sort();
            // calculate similarity scores
            for (int i = 0; i < pdbIdList.Count; i++)
            {
                InterfaceChains[] crystInterfacesI = entryCrystInterfacesDict[pdbIdList[i]];
                for (int j = i + 1; j < pdbIdList.Count; j++)
                {
                    InterfaceChains[] crystInterfacesJ = entryCrystInterfacesDict[pdbIdList[j]];
                    InterfacePairInfo[] compInfo = interfaceComp.CompareInterfacesBetweenCrystals(crystInterfacesI, crystInterfacesJ);
                    foreach (InterfacePairInfo pairInfo in compInfo)
                    {
                        DataRow dataRow = difEntryInterfaceCompTable.NewRow();
                        dataRow["PdbID1"] = pairInfo.interfaceInfo1.pdbId;
                        dataRow["InterfaceID1"] = pairInfo.interfaceInfo1.interfaceId;
                        dataRow["PdbID2"] = pairInfo.interfaceInfo2.pdbId;
                        dataRow["InterfaceID2"] = pairInfo.interfaceInfo2.interfaceId;
                        dataRow["Qscore"] = pairInfo.qScore;
                        difEntryInterfaceCompTable.Rows.Add(dataRow);
                    }
                }
            }

            return difEntryInterfaceCompTable;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="entryCrystInterfacesDict"></param>
        /// <returns></returns>
        private DataTable CompareDiffEntryInterfaces(Dictionary<string, InterfaceChains[]> entryCrystInterfacesDict, DataTable seqAlignmentTable)
        {
            DataTable difEntryInterfaceCompTable = InitializeDifEntryInterfaceCompTable();
            List<string> pdbIdList = new List<string>(entryCrystInterfacesDict.Keys);
            pdbIdList.Sort();
            UpdateResidueNumbersByInputAlignments(entryCrystInterfacesDict, seqAlignmentTable);
            // calculate similarity scores
            for (int i = 0; i < pdbIdList.Count; i++)
            {
                InterfaceChains[] crystInterfacesI = entryCrystInterfacesDict[pdbIdList[i]];
                for (int j = i + 1; j < pdbIdList.Count; j++)
                {
                    InterfaceChains[] crystInterfacesJ = entryCrystInterfacesDict[pdbIdList[j]];
                    InterfacePairInfo[] compInfo = interfaceComp.CompareInterfacesBetweenCrystals(crystInterfacesI, crystInterfacesJ);
                    foreach (InterfacePairInfo pairInfo in compInfo)
                    {
                        DataRow dataRow = difEntryInterfaceCompTable.NewRow();
                        dataRow["PdbID1"] = pairInfo.interfaceInfo1.pdbId;
                        dataRow["InterfaceID1"] = pairInfo.interfaceInfo1.interfaceId;
                        dataRow["PdbID2"] = pairInfo.interfaceInfo2.pdbId;
                        dataRow["InterfaceID2"] = pairInfo.interfaceInfo2.interfaceId;
                        dataRow["Qscore"] = pairInfo.qScore;
                        difEntryInterfaceCompTable.Rows.Add(dataRow);
                    }
                }
            }

            return difEntryInterfaceCompTable;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="entryCrystInterfacesDict"></param>
        /// <param name="seqAlignmentTable"></param>
        private void UpdateResidueNumbersByInputAlignments (Dictionary<string, InterfaceChains[]> entryCrystInterfacesDict, DataTable seqAlignmentTable)
        {
            Dictionary<string, Dictionary<string, string>> entryEntityResidueMapDict = new Dictionary<string, Dictionary<string, string>>();
            foreach (string pdbId in entryCrystInterfacesDict.Keys)
            {
                for (int i = 0; i < entryCrystInterfacesDict[pdbId].Length; i ++)
                {
                    UpdateCrystInterfaceResidueNumbers(entryCrystInterfacesDict[pdbId][i], seqAlignmentTable, entryEntityResidueMapDict);
                }
            }
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="crystInterface"></param>
        /// <param name="seqAlignmentTable"></param>
        /// <param name="longestSeqAlignment"></param>
        /// <param name="entryEntityResidueMapDict"></param>
        private void UpdateCrystInterfaceResidueNumbers (InterfaceChains crystInterface, DataTable seqAlignmentTable, Dictionary<string, Dictionary<string, string>> entryEntityResidueMapDict)
        {
            Dictionary<string, string> residueNumMapDict1 = GetResidueAlignMapDict (crystInterface.pdbId, crystInterface.entityId1.ToString (), seqAlignmentTable, entryEntityResidueMapDict);
            UpdateAtomResidueNumbers(crystInterface.chain1, residueNumMapDict1);

            Dictionary<string, string> residueNumMapDict2 = GetResidueAlignMapDict(crystInterface.pdbId, crystInterface.entityId2.ToString(), seqAlignmentTable, entryEntityResidueMapDict);
            UpdateAtomResidueNumbers(crystInterface.chain2, residueNumMapDict2);

            crystInterface.GetInterfaceResidueDist();
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="atoms"></param>
        /// <param name="residueNumberMapDict"></param>
        private void UpdateAtomResidueNumbers (AtomInfo[] atoms, Dictionary<string, string> residueNumberMapDict)
        {
            for (int i = 0; i < atoms.Length; i ++)
            {
                if (residueNumberMapDict.ContainsKey (atoms[i].seqId ))
                {
                    atoms[i].seqId = residueNumberMapDict[atoms[i].seqId];
                }
            }
        }    

        /// <summary>
        /// map residues to alignment positions
        /// </summary>
        /// <param name="pdbId"></param>
        /// <param name="entityId"></param>
        /// <param name="seqAlignmentTable"></param>
        /// <param name="longestSeqAlignment"></param>
        /// <param name="entryEntityResidueMapDict"></param>
        /// <returns></returns>
        private Dictionary<string, string> GetResidueAlignMapDict (string pdbId, string entityId, DataTable seqAlignmentTable, Dictionary<string, Dictionary<string, string>> entryEntityResidueMapDict)
        {
            if (entryEntityResidueMapDict.ContainsKey (pdbId + entityId))
            {
                return entryEntityResidueMapDict[pdbId + entityId];
            }
            DataRow[] seqAlignRows = seqAlignmentTable.Select(string.Format ("PdbID = '{0}' AND EntityID = '{1}'", pdbId, entityId));
            string seqAlignment = "";
            Dictionary<string, string> entityResidueMapDict = new Dictionary<string, string>();
            int resNo = 0;
            if (seqAlignRows.Length > 0)
            {
                seqAlignment = seqAlignRows[0]["SeqAlignment"].ToString();
                for (int i = 0; i < seqAlignment.Length; i++)
                {                   
                    if (seqAlignment[i] != '-')
                    {
                        resNo++;
                        entityResidueMapDict.Add(resNo.ToString(), (i+1).ToString());
                    }                  
                }
            }
            entryEntityResidueMapDict.Add(pdbId + entityId, entityResidueMapDict);
            return entityResidueMapDict;
        }
        #endregion

        #region clustering
        /// <summary>
        /// 
        /// </summary>
        /// <param name="difEntryInterfaceCompTable"></param>
        /// <param name="sameEntryInterfaceCompTable"></param>
        public DataTable ClusterInterfaces(DataTable difEntryInterfaceCompTable, DataTable sameEntryInterfaceCompTable, Dictionary<string, string> entrySgAsuHash, Dictionary<string, string[]> entryRepHomoInterfacesDict)
        {
            Dictionary<string, int> interfaceIndexHash = interfaceCluster.GetInterfacesIndexes(difEntryInterfaceCompTable);
          
            interfaceCluster.GetPdbInterfacesFromIndex(interfaceIndexHash);

            bool canClustered = true;
            double[,] distMatrix = interfaceCluster.CreateDistMatrix(sameEntryInterfaceCompTable, difEntryInterfaceCompTable,
                interfaceIndexHash, out canClustered);
//            Dictionary<int, int[]> clusterHash = interfaceCluster.ClusterThisGroupInterfaces(distMatrix, interfaceIndexHash);
            Dictionary<int, int[]> clusterHash = interfaceCluster.ClusterThisGroupInterfaces(distMatrix);

            DataTable clusterTable = null;
            if (clusterHash.Count > 0)
            {
                clusterTable = AssignDataToTable(clusterHash, interfaceIndexHash, entrySgAsuHash, entryRepHomoInterfacesDict);
            }

            return clusterTable;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="interfaceIndexHash"></param>
        public string[] GetPdbInterfacesFromIndex(Dictionary<string, int> interfaceIndexHash)
        {
           string[] pdbInterfaces = new string[interfaceIndexHash.Count];
            foreach (string pdbInterface in interfaceIndexHash.Keys)
            {
                pdbInterfaces[(int)interfaceIndexHash[pdbInterface]] = pdbInterface;
            }

            return pdbInterfaces;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="pdbIds"></param>
        /// <param name="asuTable"></param>
        /// <param name="entryInfoTable"></param>
        /// <returns></returns>
        private Dictionary<string, string> GetEntrySpacegroupAsuDict (string[] pdbIds, DataTable asuTable, DataTable entryInfoTable)
        {
            Dictionary<string, string> entrySgAsuDict = new Dictionary<string, string>();
            string asu = "";
            string spaceGroup = "";
            
            foreach (string pdbId in pdbIds)
            {
                asu = GetEntryAsuFormat(pdbId, asuTable);
                spaceGroup = GetEntrySpaceGroup(pdbId, entryInfoTable);
                entrySgAsuDict.Add(pdbId, spaceGroup + "_" + asu);
            }
            return entrySgAsuDict;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="entrySgAsuDict"></param>
        /// <returns></returns>
        private Dictionary<string, int> GetSgAsuCfIdDict (Dictionary<string, string> entrySgAsuDict)
        {
            Dictionary<string, int> sgAsuCfIdDict = new Dictionary<string, int>();
            List<string> sgAsuList = new List<string>();
            foreach (string pdbId in entrySgAsuDict.Keys)
            {
               if (! sgAsuList.Contains (entrySgAsuDict[pdbId]))
               {
                   sgAsuList.Add(entrySgAsuDict[pdbId]);
               }
            }
            sgAsuList.Sort();
             int cfId = 1;
            foreach (string sgAsu in sgAsuList)
            {
                sgAsuCfIdDict.Add(sgAsu, cfId);
                cfId++;
            }
            return sgAsuCfIdDict;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="pdbId"></param>
        /// <param name="asuTable"></param>
        /// <returns></returns>
        private string GetEntryAsuFormat (string pdbId, DataTable asuTable)
        {
            DataRow[] asuRows = asuTable.Select(string.Format ("PdbID = '{0}' AND PolymerType = 'polypeptide'", pdbId));
            Dictionary<string, int> protEntityCountDict = new Dictionary<string, int>();
            string entityId = "";
            foreach (DataRow asuRow in asuRows)
            {
                entityId = asuRow["EntityID"].ToString ();
                if (protEntityCountDict.ContainsKey (entityId))
                {
                    protEntityCountDict[entityId] = protEntityCountDict[entityId] + 1;
                }
                else
                {
                    protEntityCountDict.Add(entityId, 1);
                }
            }
            string asuFormat = buQuery.GetAbcFormatFromEntityHash(protEntityCountDict);
            return asuFormat;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="pdbId"></param>
        /// <param name="entryInfoTable"></param>
        /// <returns></returns>
        private string GetEntrySpaceGroup(string pdbId, DataTable entryInfoTable)
        {
            DataRow[] entryRows = entryInfoTable.Select(string.Format("PdbID = '{0}'", pdbId));
            if (entryRows.Length > 0)
            {
                return entryRows[0]["SpaceGroup"].ToString();
            }
            return "-";
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="clusterHash"></param>
        /// <param name="interfaceIndexHash"></param>
        /// <returns></returns>
        private DataTable AssignDataToTable(Dictionary<int, int[]> clusterHash, Dictionary<string, int> interfaceIndexHash, Dictionary<string, string> entrySgAsuDict,
            Dictionary<string, string[]> entryRepHomoInterfacesDict)
        {
            DataTable clusterTable = InitializeClusterTable();
            int groupId = 111111;
            string[] pdbInterfaces = GetPdbInterfacesFromIndex(interfaceIndexHash);
            Dictionary<string, int> sgAsuCfIdDict = GetSgAsuCfIdDict(entrySgAsuDict);
            int[][] clusterNumEntryLists = SortClustersByNumEntries(clusterHash, interfaceIndexHash);
            int[] clusterIds = clusterNumEntryLists[0];
            int[] clusterEntryNumbers = clusterNumEntryLists[1]; 
            string[] sgAsuFields = null;
            string sgAsu = "";
            List<string> clusterHomoInterfaceList = new List<string>();
            string[] pdbInterfaceFields = null;
            
            int realClusterId = 1;
            for (int i = 0; i < clusterIds.Length; i++)
            {
                if (clusterEntryNumbers[i] <= 1)
                {
                    continue;
                }
                int[] interfaceList = clusterHash[clusterIds[i]];
                clusterHomoInterfaceList.Clear();
                foreach (int interfaceIndex in interfaceList)
                {
                    string pdbInterfaceString = pdbInterfaces[interfaceIndex];
                    pdbInterfaceFields = pdbInterfaceString.Split('_');
                    sgAsu = entrySgAsuDict[pdbInterfaceFields[0]];
                    DataRow clusterRow = clusterTable.NewRow();
                    clusterRow["GroupSeqID"] = groupId;
                    clusterRow["CfGroupID"] = sgAsuCfIdDict[sgAsu];
                    clusterRow["ClusterID"] = realClusterId;
                    clusterRow["PdbID"] = pdbInterfaceFields[0];
                    clusterRow["InterfaceID"] = pdbInterfaceFields[1];
                    sgAsuFields = entrySgAsuDict[pdbInterfaceFields[0]].Split('_');
                    clusterRow["SpaceGroup"] = sgAsuFields[0];
                    clusterRow["ASU"] = sgAsuFields[1];
                    clusterTable.Rows.Add(clusterRow);

                    if (entryRepHomoInterfacesDict.ContainsKey(pdbInterfaceString))
                    {

                        foreach (string homoInterface in entryRepHomoInterfacesDict[pdbInterfaceString])
                        {
                            if (!clusterHomoInterfaceList.Contains(homoInterface))
                            {
                                clusterHomoInterfaceList.Add(homoInterface);
                            }
                        }
                    }
                }

                foreach (string homoInterface in clusterHomoInterfaceList)
                {
                    pdbInterfaceFields = homoInterface.Split('_');
                    sgAsu = entrySgAsuDict[pdbInterfaceFields[0]];
                    DataRow homoClusterRow = clusterTable.NewRow();
                    homoClusterRow["GroupSeqID"] = groupId;
                    homoClusterRow["CfGroupID"] = sgAsuCfIdDict[sgAsu];
                    homoClusterRow["ClusterID"] = realClusterId;
                    homoClusterRow["PdbID"] = pdbInterfaceFields[0];
                    homoClusterRow["InterfaceID"] = pdbInterfaceFields[1];
                    sgAsuFields = entrySgAsuDict[pdbInterfaceFields[0]].ToString().Split('_');
                    homoClusterRow["SpaceGroup"] = sgAsuFields[0];
                    homoClusterRow["ASU"] = sgAsuFields[1];
                    clusterTable.Rows.Add(homoClusterRow);
                }
                realClusterId++;
            }
            return clusterTable;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="clusterHash"></param>
        /// <param name="interfaceIndexHash"></param>
        /// <returns></returns>
        private int[][] SortClustersByNumEntries (Dictionary<int, int[]> clusterHash, Dictionary<string, int> interfaceIndexHash)
        {
            List<string> entryList = new List<string> ();
            string[] pdbInterfaces = GetPdbInterfacesFromIndex(interfaceIndexHash);
            List<int> clusterIdList = new List<int>(clusterHash.Keys);
            int[] clusterIds = clusterIdList.ToArray();
            int[] clusterEntryNumbers = new int[clusterIds.Length];
            for (int i = 0; i < clusterIds.Length; i ++ )
            {
                int[] interfaceList = clusterHash[clusterIds[i]];

                entryList.Clear();
                foreach (int interfaceIndex in interfaceList)
                {
                    string pdbInterfaceString = pdbInterfaces[interfaceIndex];
                    string[] pdbInterfaceFields = pdbInterfaceString.Split('_');
                    if (!entryList.Contains(pdbInterfaceFields[0]))
                    {
                        entryList.Add(pdbInterfaceFields[0]);
                    }
                }
                clusterEntryNumbers[i] = entryList.Count;
            }
            for (int i = 0; i < clusterIds.Length; i ++)
            {
                for (int j = i + 1; j < clusterIds.Length; j ++)
                {
                    if (clusterEntryNumbers[i] < clusterEntryNumbers[j])
                    {
                        int temp = clusterEntryNumbers[i];
                        clusterEntryNumbers[i] = clusterEntryNumbers[j];
                        clusterEntryNumbers[j] = temp;

                        int clusterIdTemp = clusterIds[i];
                        clusterIds[i] = clusterIds[j];
                        clusterIds[j] = clusterIdTemp;
                    }
                }
            }
            int[][] clusterIdNumEntryLists = new int[2][];
            clusterIdNumEntryLists[0] = clusterIds;
            clusterIdNumEntryLists[1] = clusterEntryNumbers.ToArray();

            return clusterIdNumEntryLists;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        private DataTable InitializeClusterTable ()
        {
            DataTable clusterTable = new DataTable("Clusters");
            string[] clusterCols = {"GroupSeqID", "CfGroupID", "ClusterID", "PdbID", "InterfaceID", "SpaceGroup", "ASU"};
            foreach (string col in clusterCols)
            {
                clusterTable.Columns.Add(new DataColumn(col));
            }
            return clusterTable;
        }
        #endregion

        #region compile cluster interfaces including pymol scripts
        /// <summary>
        /// 
        /// </summary>
        /// <param name="userGroupName"></param>
        /// <param name="clusterTable"></param>
        /// <param name="clusterFileDir"></param>
        /// <param name="interfaceFileDir"></param>
        private void CompileClusterCoordinates (string userGroupName, DataTable clusterTable, string clusterFileDir, string interfaceFileDir)
        {           
            string pmlScriptFileName = userGroupName;
            
            Dictionary<int, List<string>> clusterInterfacesDict = new Dictionary<int, List<string>>();
            int clusterId = 0;
            string clusterInterface = "";
            foreach (DataRow interfaceRow in clusterTable.Rows)
            {
                clusterId = Convert.ToInt32(interfaceRow["ClusterID"].ToString ());
                clusterInterface = interfaceRow["PdbID"].ToString() + "_" + interfaceRow["InterfaceID"].ToString();
                if (clusterInterfacesDict.ContainsKey (clusterId))
                {
                    clusterInterfacesDict[clusterId].Add(clusterInterface); 
                }
                else
                {
                    List<string> clusterInterfaceList = new List<string>();
                    clusterInterfaceList.Add(clusterInterface);
                    clusterInterfacesDict.Add(clusterId, clusterInterfaceList);
                }
            }
            List<int> clusterIdList = new List<int>(clusterInterfacesDict.Keys);
            clusterIdList.Sort();
            string clusterFile = "";
              List<string> clusterFileList = new List<string> ();
            foreach (int lsClusterId in clusterIdList)
            {
                try
                {
                    clusterFile = CompressGroupClusterInterfaceFiles(lsClusterId, userGroupName, clusterTable, interfaceFileDir, clusterFileDir);
                    clusterFileList.Add(clusterFile);
                }
                catch (Exception ex)
                {
                    Console.WriteLine("Compress cluster coordinates warning: " + lsClusterId.ToString () + " " + ex.Message);
                }
            }

            // tar cluster files to group
            string groupTarFileName = userGroupName + ".tar";
            fileCompress.RunTar(groupTarFileName, clusterFileList.ToArray(), clusterFileDir, false);
        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="superGroupId"></param>
        /// <param name="clusterId"></param>
        /// <param name="clusterReverseFilesHash"></param>
        /// <returns></returns>
        public string CompressGroupClusterInterfaceFiles(int clusterId, string groupName, DataTable clusterTable, string interfaceFileDir, string clusterFileDir)
        {
            DataRow[] clusterInterfaceRows = clusterTable.Select(string.Format ("ClusterID = '{0}'", clusterId));
            List<string> entryList = new List<string>();
            List<string> entryInterfaceList = new List<string>();
            foreach (DataRow entryRow in clusterInterfaceRows)
            {
                string pdbId = entryRow["PdbID"].ToString();
                /*          if (!entryList.Contains(pdbId))
                          {
                              entryList.Add(pdbId);*/
                entryInterfaceList.Add(pdbId + "_" + entryRow["InterfaceID"].ToString());
                //         }
            }
            string interfaceFile = "";
            List<string> clusterInterfaceFileList = new List<string>();
            foreach (string entryInterface in entryInterfaceList)
            {
                interfaceFile = Path.Combine(interfaceFileDir, entryInterface + ".cryst.gz");
                try
                {
                    File.Copy(interfaceFile, Path.Combine(clusterFileDir, entryInterface + ".cryst.gz"), true);
                }
                catch (Exception ex)
                {
                    Console.WriteLine (ex.Message);
                    continue;
                }
                ParseHelper.UnZipFile(Path.Combine(clusterFileDir, entryInterface + ".cryst.gz"));
                clusterInterfaceFileList.Add(entryInterface + ".cryst");
            }

            string pmlScriptFileName = groupName + "_" + clusterId;
            string[] pymolScriptFiles = pmlScriptWriter.FormatChainInterfacePymolScriptFiles(pmlScriptFileName, clusterFileDir, clusterInterfaceFileList.ToArray ());

            // add pymol script files into cluster tar filed
            clusterInterfaceFileList.AddRange(pymolScriptFiles);

            string clusterFileName = groupName + "_" + clusterId + ".tar.gz";
            string clusterFile = fileCompress.RunTar(clusterFileName, clusterInterfaceFileList.ToArray(), clusterFileDir, true);

            foreach (string clusterInterfaceFile in clusterInterfaceFileList)
            {
                File.Delete(Path.Combine(clusterFileDir, clusterInterfaceFile));
            }
            return clusterFile;
        }
        #endregion

        #region print to text files
        /// <summary>
        /// 
        /// </summary>
        /// <param name="asuTables"></param>
        /// <param name="textDataDir"></param>
        private void PrintEntryAsuInfoToTextFiles (DataTable[] asuTables, string textDataDir)
        {
            string asuFile = "";
            string headerLine = "";
            for (int i = 0; i < asuTables.Length; i ++)
            {
                asuFile = Path.Combine(textDataDir, asuTables[i].TableName + ".txt");
                StreamWriter dataWriter = new StreamWriter(asuFile);
                headerLine = "";
                foreach (DataColumn dCol in asuTables[i].Columns)
                {
                    headerLine += dCol.ColumnName + "\t";
                }
                dataWriter.WriteLine(headerLine.TrimEnd('\t'));
                dataWriter.WriteLine();
                dataWriter.WriteLine(ParseHelper.FormatDataRows(asuTables[i].Select()));
                dataWriter.Close();
            }
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="entryCrystInterfacesDict"></param>
        /// <param name="interfaceTextFile"></param>
        private void PrintEntryCrystInterfacesToFile (Dictionary<string, InterfaceChains[]> entryCrystInterfacesDict, string interfaceTextFile, DataTable asuTable)
        {
            StreamWriter dataWriter = new StreamWriter(interfaceTextFile);
            dataWriter.WriteLine("PdbID\tInterfaceID\tAsymID1\tEntityID1\tAuthorChain1\tSymmetryString1\tAsymID2\tEntityID2\tAuthorChain2\tSymmetryString2\tSurfaceArea");
            string dataLine = "";
            List<string> entryList = new List<string>(entryCrystInterfacesDict.Keys);
            entryList.Sort();
            foreach (string pdbId in entryList)
            {
                foreach (InterfaceChains crystInterface in entryCrystInterfacesDict[pdbId])
                {
                    dataLine = FormatInterfaceLine (crystInterface, asuTable);
                    dataWriter.WriteLine(dataLine);
                }
            }
            dataWriter.Close();
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="crystInterface"></param>
        /// <returns></returns>
        private string FormatInterfaceLine(InterfaceChains crystInterface, DataTable asuTable)
        {
            int entityId1 = -1;
            string authorChain1 = "";
            string symOpStr1 = crystInterface.firstSymOpString;
            string asymChain1 = symOpStr1.Substring(0, symOpStr1.IndexOf("_"));
            string symString1 = symOpStr1.Substring(symOpStr1.IndexOf("_") + 1, symOpStr1.Length - symOpStr1.IndexOf("_") - 1);
            GetAuthEntityInfoForAsymChain(asymChain1, asuTable, out entityId1, out authorChain1);
            int entityId2 = -1;
            string authorChain2 = "";
            string symOpStr2 = crystInterface.secondSymOpString;
            string asymChain2 = symOpStr2.Substring(0, symOpStr2.IndexOf("_"));
            string symString2 = symOpStr2.Substring(symOpStr2.IndexOf("_") + 1, symOpStr2.Length - symOpStr2.IndexOf("_") - 1);
            if (asymChain1 == asymChain2)
            {
                entityId2 = entityId1;
                authorChain2 = authorChain1;
            }
            else
            {
                GetAuthEntityInfoForAsymChain(asymChain2, asuTable, out entityId2, out authorChain2);
            }

            string dataLine = crystInterface.pdbId + "\t" + crystInterface.interfaceId + "\t" +
                asymChain1 + "\t" + entityId1 + "\t" + authorChain1 + "\t" + symString1 + "\t" +
                asymChain2 + "\t" + entityId2 + "\t" + authorChain2 + "\t" + symString2 + "\t" +
                crystInterface.surfaceArea;

            return dataLine;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="sameEntryInterfaceCompTable"></param>
        /// <param name="sameEntryInterfaceCompFile"></param>
        private void PrintSameEntryInterfaceCompToTextFile (DataTable sameEntryInterfaceCompTable, string sameEntryInterfaceCompFile)
        {
            StreamWriter dataWriter = new StreamWriter(sameEntryInterfaceCompFile);
            dataWriter.WriteLine("PdbID\tInterfaceID1\tInterfaceID2\tQscore");
            foreach (DataRow compRow in sameEntryInterfaceCompTable.Rows)
            {
                dataWriter.WriteLine(compRow["PdbID"] + "\t" + compRow["InterfaceID1"] + "\t" +
                    compRow["InterfaceID2"] + "\t" + compRow["Qscore"]);
            }
            dataWriter.Close();
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="sameEntryInterfaceCompTable"></param>
        /// <param name="sameEntryInterfaceCompFile"></param>
        private void PrintDiffEntryInterfaceCompToTextFile(DataTable diffEntryInterfaceCompTable, string diffEntryInterfaceCompFile)
        {
            StreamWriter dataWriter = new StreamWriter(diffEntryInterfaceCompFile);
            dataWriter.WriteLine("PdbID1\tInterfaceID1\tPdbID2\tInterfaceID2\tQscore");
            foreach (DataRow compRow in diffEntryInterfaceCompTable.Rows)
            {
                dataWriter.WriteLine(compRow["PdbID1"] + "\t" + compRow["InterfaceID1"] + "\t" +
                    compRow["PdbID2"] + "\t" + compRow["InterfaceID2"] + "\t" + compRow["Qscore"]);
            }
            dataWriter.Close();
        }

        /// <summary>
        ///  
        /// </summary>
        /// <param name="clusterTable"></param>
        /// <param name="clusterTextFile"></param>
        private void PrintClusterTableToTextFile (DataTable clusterTable, string clusterTextFile)
        {
            StreamWriter dataWriter = new StreamWriter(clusterTextFile);
            dataWriter.WriteLine ("ClusterID\tCF ID\tPdbID\tInterfaceID\tSpaceGroup\tASU");
            
            foreach (DataRow interfaceRow in clusterTable.Rows)
            {
                dataWriter.WriteLine(interfaceRow["ClusterID"] + "\t" + interfaceRow["CfGroupID"] + "\t" + 
                    interfaceRow["PdbID"] + "\t" + interfaceRow["InterfaceID"]+ "\t" + 
                    interfaceRow["SpaceGroup"] + "\t" + interfaceRow["ASU"]);
            }
            dataWriter.Close();
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="asymChain"></param>
        /// <param name="chainInfoTable"></param>
        /// <param name="entityId"></param>
        /// <param name="authorChain"></param>
        private void GetAuthEntityInfoForAsymChain(string asymChain, DataTable chainInfoTable, out int entityId, out string authorChain)
        {
            entityId = -1;
            authorChain = "-";
            DataRow[] chainRows = chainInfoTable.Select(string.Format("AsymID = '{0}'", asymChain));
            if (chainRows.Length > 0)
            {
                entityId = Convert.ToInt32(chainRows[0]["EntityID"].ToString());
                authorChain = chainRows[0]["AuthorChain"].ToString().TrimEnd();
            }          
        }
        #endregion

        #region helper functions
        /// <summary>
        /// 
        /// </summary>
        /// <param name="entryLsFile"></param>
        /// <returns></returns>
        private string[] ReadEntries(string entryLsFile)
        {
            List<string> entryList = new List<string>();
            StreamReader entryReader = new StreamReader(entryLsFile);
            string line = "";
            while ((line = entryReader.ReadLine()) != null)
            {
                entryList.Add(line);
            }
            entryReader.Close();
            entryList.Sort();
            return entryList.ToArray();
        }

        /// <summary>
        /// clustal omega format
        /// CLUSTAL O(1.2.4) multiple sequence alignment
        /*
        1aaaA      --------------------------------MSE----NPLLERARRFLSALRHCQVLG	24
        1bbbB      ---------------------------------MS----STALEMASRFVNRSPFNRWLG	23
                                                                                  

        1aaaA      --L-TVEAADEK----------GLTLRLPYSQAIIGNPESGVVHGGAITTLMDTTCGIST	71
        1bbbB      --M-SVLEAGEQ----------GIVLGIKWREELISSPEIRSTHGGILATLVDAAGDYAV	70 
        */
        /// OR
        /*
         * 1aaa1  --------------------------------MSE----NPLLERARRFLSALRHCQVLG--L-TVEAADEK----------GLTLRLPYSQAIIGNPESGVVHGGAITTLMDTTCGIST
         * 1bbb2  ---------------------------------MS----STALEMASRFVNRSPFNRWLG--M-SVLEAGEQ----------GIVLGIKWREELISSPEIRSTHGGILATLVDAAGDYAV
         * */
        /// sequence name must be PdbID + entity ID (integer) or chain ID 
        /// AND residue numbers must be in 1:N
        /// must be separated by white space between sequence name and sequence alignment
        /// gap must be '-'
        /// </summary>
        /// <param name="alignFile"></param>
        /// <returns></returns>
        private DataTable ReadAlignmentFile(string alignFile, string[] pdbIds)
        {
            DataTable seqAlignmentTable = InitializeAlignmentTable();
            string line = "";
            StreamReader dataReader = new StreamReader  (alignFile);
            string pdbId = "";
            int entityId = 0;
            string chainId = "";
            string entityOrChainField = "";
            Dictionary<string, string> entryAlignmentDict = new Dictionary<string, string>();
            while ((line = dataReader.ReadLine ()) != null)
            {
                string[] lineFields = ParseHelper.SplitPlus(line, ' ');
                if (lineFields.Length >= 2 && lineFields[0].Length > 4)
                {
                    pdbId = lineFields[0].Substring(0, 4);
                    if (Array.IndexOf(pdbIds, pdbId) > -1)
                    {                        
                        string[] alignFields = lineFields[1].Split('\t'); // alignment: alignFields[0]
                        if (entryAlignmentDict.ContainsKey (lineFields[0]))
                        {
                            entryAlignmentDict[lineFields[0]] = entryAlignmentDict[lineFields[0]] + alignFields[0];
                        }
                        else
                        {
                            entryAlignmentDict.Add(lineFields[0], alignFields[0]);
                        }
                    }
                }
            }
            dataReader.Close();
           
            foreach (string entryName in entryAlignmentDict.Keys)
            {
                pdbId = entryName.Substring(0, 4);
                chainId = "-";
                entityId = -1;
                entityOrChainField = entryName.Substring(4, entryName.Length - 4);
                if (!Int32.TryParse(entityOrChainField, out entityId))
                {
                    chainId = entityOrChainField;
                    entityId = -1;
                }
                DataRow dataRow = seqAlignmentTable.NewRow();
                dataRow["PdbID"] = pdbId;
                dataRow["EntityID"] = entityId;
                dataRow["ChainID"] = chainId;
                dataRow["SeqAlignment"] = entryAlignmentDict[entryName];
                seqAlignmentTable.Rows.Add(dataRow);
            }
            return seqAlignmentTable;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        private DataTable InitializeAlignmentTable ()
        {
            DataTable alignmentTable = new DataTable("SeqAlignments");
            string[] alignColumns = {"PdbID", "EntityID", "ChainID", "SeqAlignment"};
            foreach (string col in alignColumns)
            {
                alignmentTable.Columns.Add(new DataColumn(col));
            }
            return alignmentTable;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="seqAlignmentTable"></param>
        /// <param name="asuTable"></param>
        private void FillAlignmentTableEntityChainIds (DataTable seqAlignmentTable, DataTable asuTable)
        {
            string chainId = "";
            string entityId = "";
            for (int i = 0; i < seqAlignmentTable.Rows.Count; i++)
            {
                if (seqAlignmentTable.Rows[i]["EntityID"].ToString () == "-1")
                {
                    entityId = GetEntityIdFromAuthChain(seqAlignmentTable.Rows[i]["PdbID"].ToString(), seqAlignmentTable.Rows[i]["ChainID"].ToString(), asuTable);
                    seqAlignmentTable.Rows[i]["EntityID"] = entityId;
                }
                else
                {
                    chainId = GetAuthIdFromEntityId(seqAlignmentTable.Rows[i]["PdbID"].ToString(), seqAlignmentTable.Rows[i]["EntityID"].ToString(), asuTable);
                    seqAlignmentTable.Rows[i]["ChainID"] = chainId;
                }
            }
            seqAlignmentTable.AcceptChanges();
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="pdbId"></param>
        /// <param name="chainId"></param>
        /// <param name="asuTable"></param>
        /// <returns></returns>
        private string GetEntityIdFromAuthChain (string pdbId, string chainId, DataTable asuTable)
        {
            DataRow[] dataRows = asuTable.Select(string.Format("PdbID = '{0}' AND AuthorChain = '{1}'", pdbId, chainId));
            if (dataRows.Length > 0)
            {
                return dataRows[0]["EntityID"].ToString();
            }
            return "-1";
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="pdbId"></param>
        /// <param name="chainId"></param>
        /// <param name="asuTable"></param>
        /// <returns></returns>
        private string GetEntityIdFromAsymId (string pdbId, string chainId, DataTable asuTable)
        {
            DataRow[] dataRows = asuTable.Select(string.Format("PdbID = '{0}' AND AsymID = '{1}'", pdbId, chainId));
            if (dataRows.Length > 0)
            {
                return dataRows[0]["EntityID"].ToString();
            }
            return "-1";
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="pdbId"></param>
        /// <param name="entityId"></param>
        /// <param name="asuTable"></param>
        /// <returns></returns>
        private string GetAuthIdFromEntityId (string pdbId, string entityId, DataTable asuTable)
        {
            DataRow[] dataRows = asuTable.Select(string.Format("PdbID = '{0}' AND EntityID = '{1}'", pdbId, entityId));
            string authId = "";
            if (dataRows.Length > 0)
            {
                authId = dataRows[0]["AuthorChain"].ToString();
            }
            return authId;
        }
        #endregion
    }
}
