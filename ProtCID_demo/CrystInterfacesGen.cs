using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Data;
using CrystalInterfaceLib.FileParser;
using CrystalInterfaceLib.Contacts;
using CrystalInterfaceLib.Crystal;
using CrystalInterfaceLib.BuIO;
using CrystalInterfaceLib.ProtInterfaces;
using InterfaceClusterLib.InterfaceProcess;

namespace ProtCID_demo
{
    public class CrystInterfacesGen
    {
        #region member variables
        private InterfaceFileWriter interfaceWriter = new InterfaceFileWriter();
        public InterfaceWriter interfaceChainWriter = new InterfaceWriter();
        private InterfaceAsa interfaceAsa = new InterfaceAsa();
 //       private StreamWriter logWriter = new StreamWriter("protcid_demo_log.txt");
        #endregion

        #region interfaces from crystal
        /// <summary>
        /// 
        /// </summary>
        /// <param name="pdbId"></param>
        /// <param name="crystFileDir"></param>
        public InterfaceChains[] GenerateInterfacesFromCryst(string pdbId, string crystFileDir, string interfaceFileDir, DataTable asuTable, string residueNumberingType)
        {
            CrystalXmlBuilder xmlBuilder = new CrystalXmlBuilder();
            string pdbXmlFile = Path.Combine(crystFileDir, pdbId + ".xml");
            string coordXmlFile = Path.Combine(crystFileDir, pdbId + "_coord.xml"); ;
            if (!File.Exists(coordXmlFile))
            {
                coordXmlFile = xmlBuilder.WritePdbXmlToMyCoordXmlFile(pdbXmlFile);
            }
            InterfaceChains[] entryInterfaces = GenerateInterfacesFromCryst(coordXmlFile, interfaceFileDir, asuTable, residueNumberingType);
            if (entryInterfaces == null)
            {
                Console.WriteLine(pdbId + " no interaces in crystal");
                return null;
            }
            return entryInterfaces;       
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="coordXmlFile"></param>
        public InterfaceChains[] GenerateInterfacesFromCryst(string coordXmlFile, string interfaceFileDir, DataTable asuTable, string residueNumberingType)
        {
            ContactInCrystal contactInCrystal = new ContactInCrystal();
            int numOfChainsInUnitCell = contactInCrystal.FindInteractChains(coordXmlFile);
            InterfaceChains[] entryInterfaces = contactInCrystal.InterfaceChainsList;
            if (entryInterfaces == null)
            {
                return null;
            }
            if (residueNumberingType == "author" || residueNumberingType == "auth")
            {
                UpdateEntryInterfaceResidueNumberingToAuth(entryInterfaces);
            }

            List<ProtInterfaceInfo> interfaceList = new List<ProtInterfaceInfo>();
            FileInfo fileInfo = new FileInfo(coordXmlFile);
            string pdbId = fileInfo.Name.Substring(0, 4);
            CalculateInterfaceSurfacrAreas(pdbId, entryInterfaces);

            Dictionary<string, int> asymChainEntityHash = GetAsymChainEntityHash(pdbId, asuTable);
            string asymChain = "";
            foreach (InterfaceChains thisInterface in entryInterfaces)
            {
                asymChain = GetAsymChainFromSymOpString(thisInterface.firstSymOpString);
                thisInterface.entityId1 = (int)asymChainEntityHash[asymChain];
                asymChain = GetAsymChainFromSymOpString(thisInterface.secondSymOpString);
                thisInterface.entityId2 = (int)asymChainEntityHash[asymChain];

                ProtInterfaceInfo interfaceInfo = GetInterfaceInfo(thisInterface, asuTable);
                interfaceInfo.ASA = thisInterface.surfaceArea;
                interfaceList.Add(interfaceInfo);
            }
            bool needLigandInfo = true;
            bool isAsuInterface = false;
            if (interfaceList.Count > 0)
            {
                ProtInterfaceInfo[] interfacesInfos = new ProtInterfaceInfo[interfaceList.Count];
                interfaceList.CopyTo(interfacesInfos);
                try
                {
                    if (residueNumberingType == "author" || residueNumberingType == "auth")
                    {
                        interfaceWriter.GenerateInterfaceFilesInAuthNumbering
                           (pdbId, coordXmlFile, interfacesInfos, "cryst", needLigandInfo, interfaceFileDir, asuTable, isAsuInterface);
                    }
                    else
                    {
                        interfaceWriter.GenerateInterfaceFiles
                            (pdbId, coordXmlFile, interfacesInfos, "cryst", needLigandInfo, interfaceFileDir, asuTable, isAsuInterface);
                    }
                }
                catch (Exception ex)
                {
                    throw new Exception(pdbId + "generate interfaces from crystal error: " + ex.Message);
                }
            }
            return entryInterfaces;
        } 

        /// <summary>
        /// 
        /// </summary>
        /// <param name="symOpString"></param>
        /// <returns></returns>
        private string GetAsymChainFromSymOpString(string symOpString)
        {
            string[] fields = symOpString.Split('_');
            return fields[0];
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="pdbId"></param>
        /// <returns></returns>
        private Dictionary<string, int> GetAsymChainEntityHash(string pdbId, DataTable asuTable)
        {
            DataRow[] entryRows = asuTable.Select(string.Format ("PdbID = '{0}'", pdbId));
            Dictionary<string, int> asymChainEntityHash = new Dictionary<string, int>();
            int entityId = 0;
            string asymChain = "";
            foreach (DataRow chainRow in entryRows)
            {
                entityId = Convert.ToInt32(chainRow["EntityID"].ToString());
                asymChain = chainRow["AsymID"].ToString().TrimEnd();
                asymChainEntityHash.Add(asymChain, entityId);
            }
            return asymChainEntityHash;
        }
   
        /// <summary>
        /// 
        /// </summary>
        /// <param name="theInterface"></param>
        /// <returns></returns>
        private string  GetInterfaceRemark (InterfaceChains theInterface, DataTable asuTable)
        {
            string authorChain1 = "";
            string authorChain2 = "";
            string[] chainSymOpStrings1 = GetAsymChainAndSymOpString(theInterface.firstSymOpString);
            authorChain1 = GetAuthorChain(theInterface.pdbId, chainSymOpStrings1[0], asuTable);
            string[] chainSymOpStrings2 = GetAsymChainAndSymOpString(theInterface.secondSymOpString);
            if (chainSymOpStrings1[0] == chainSymOpStrings2[0])
            {
                authorChain2 = authorChain1;
            }
            else
            {
                authorChain2 = GetAuthorChain(theInterface.pdbId, chainSymOpStrings2[0], asuTable);
            }
            string remark = "Remark 300 Interface Chain A For Asymmetric Chain " +
                chainSymOpStrings1[0] + " Author Chain " + authorChain1 + " " +
                " Entity  " + theInterface.entityId1.ToString () + " " +
                " Symmetry Operator    " + chainSymOpStrings1[1] + "\r\n" +
                "Remark 300 Interface Chain B For Asymmetric Chain " +
                chainSymOpStrings2[0] + " Author Chain " + authorChain2 + " " + 
                " Entity  " + theInterface.entityId2.ToString () + " " +
                " Symmetry Operator    " + chainSymOpStrings2[1] + "\r\n" +
                "Remark 350 Interface surface area: " + theInterface.surfaceArea.ToString() + "\r\n";

            return remark;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="theInterface"></param>
        /// <returns></returns>
        private ProtInterfaceInfo GetInterfaceInfo(InterfaceChains theInterface, DataTable asuTable)
        {
            ProtInterfaceInfo interfaceInfo = new ProtInterfaceInfo();
            interfaceInfo.InterfaceId = theInterface.interfaceId;
            interfaceInfo.ASA = -1.0;
            string authorChain1 = "";
            string authorChain2 = "";
            string[] chainSymOpStrings = GetAsymChainAndSymOpString(theInterface.firstSymOpString);
            interfaceInfo.Chain1 = chainSymOpStrings[0];
            authorChain1 = GetAuthorChain(theInterface.pdbId, interfaceInfo.Chain1, asuTable);
            interfaceInfo.SymmetryString1 = chainSymOpStrings[1];
            chainSymOpStrings = GetAsymChainAndSymOpString(theInterface.secondSymOpString);
            interfaceInfo.Chain2 = chainSymOpStrings[0];
            if (interfaceInfo.Chain1 == interfaceInfo.Chain2)
            {
                authorChain2 = authorChain1;
            }
            else
            {
                authorChain2 = GetAuthorChain(theInterface.pdbId, interfaceInfo.Chain2, asuTable);
            }
            interfaceInfo.SymmetryString2 = chainSymOpStrings[1];
            interfaceInfo.Remark = "Remark 300 Interface Chain A For Asymmetric Chain " +
                interfaceInfo.Chain1 + " Author Chain " + authorChain1 + " " +
                " Entity  " + theInterface.entityId1.ToString() + " " +
                " Symmetry Operator    " + interfaceInfo.SymmetryString1 + "\r\n" +
                "Remark 300 Interface Chain B For Asymmetric Chain " +
                interfaceInfo.Chain2 + " Author Chain " + authorChain2 + " " +
                " Entity  " + theInterface.entityId2.ToString() + " " +
                " Symmetry Operator    " + interfaceInfo.SymmetryString2 + "\r\n";
            return interfaceInfo;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="asymId"></param>
        /// <param name="asuTable"></param>
        /// <returns></returns>
        public string GetAuthorChain(string pdbId, string asymId, DataTable entryInfoTable)
        {
            DataRow[] authChainRows = entryInfoTable.Select(string.Format("PdbID = '{0}' AND AsymID = '{1}'", pdbId, asymId));
            if (authChainRows.Length > 0)
            {
                return authChainRows[0]["AuthorChain"].ToString().TrimEnd();
            }
            return "-";
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="symOpString"></param>
        /// <returns></returns>
        public string[] GetAsymChainAndSymOpString(string symOpString)
        {
            string[] chainSymOpStrings = new string[2];
            string[] symOpFields = symOpString.Split('_');
            if (symOpFields.Length == 1)
            {
                chainSymOpStrings[0] = symOpFields[0];
                chainSymOpStrings[1] = "1_555";
            }
            else if (symOpFields.Length == 2)
            {
                chainSymOpStrings[0] = symOpFields[0];
                chainSymOpStrings[1] = symOpFields[1];
            }
            else if (symOpFields.Length == 3)
            {
                chainSymOpStrings[0] = symOpFields[0];
                chainSymOpStrings[1] = symOpFields[1] + "_" + symOpFields[2];
            }
            return chainSymOpStrings;
        }
        #endregion

        #region interface surface area
        /// <summary>
        /// generate interface files from the input defintion
        /// </summary>
        /// <param name="pdbId"></param>
        /// <param name="interfaceInfoList"></param>
        /// <param name="type"></param>
        /// <param name="needAsaUpdated"></param>
        /// <returns></returns>
        public void CalculateInterfaceSurfacrAreas (string pdbId, InterfaceChains[] entryInterfaces)
        {
            string type = "cryst";
            Dictionary<string, AtomInfo[]> interfaceChainsHash = new Dictionary<string,AtomInfo[]> ();
            string chainA = "";
            string chainB = "";
            double asa = 0;
            foreach (InterfaceChains thisInterface in entryInterfaces)
            {
                chainA = GetChainFromSymString(thisInterface.firstSymOpString);
                chainB = GetChainFromSymString(thisInterface.secondSymOpString);
                if (! interfaceChainsHash.ContainsKey (chainA))
                {
                    interfaceChainsHash.Add(chainA, thisInterface.chain1);
                }
                if (! interfaceChainsHash.ContainsKey (chainB))
                {
                    interfaceChainsHash.Add(chainB, thisInterface.chain2);
                }
            }
            Dictionary<string,double> chainSaHash = interfaceAsa.GetChainsSurfaceAreaInBu (pdbId, interfaceChainsHash, type);

            double chainASa = 0;
            double chainBSa = 0;
            double[] interfaceSurfaceAreas = new double[entryInterfaces.Length];
            for (int i = 0; i < entryInterfaces.Length; i ++ )
            {
                chainA = GetChainFromSymString(entryInterfaces[i].firstSymOpString);
                chainB = GetChainFromSymString(entryInterfaces[i].secondSymOpString);
                string interfaceComplexFile = interfaceChainWriter.WriteTempInterfaceToFile(pdbId, entryInterfaces[i].interfaceId, entryInterfaces[i].chain1, entryInterfaces[i].chain2, type);
                double complexAsa = interfaceAsa.ComputeInterfaceSurfaceArea(interfaceComplexFile);

                if (chainSaHash.ContainsKey(chainA))
                {
                    chainASa = (double)chainSaHash[chainA];
                }
                else
                {
                    chainASa = -1;
                }
                if (chainSaHash.ContainsKey(chainB))
                {
                    chainBSa = (double)chainSaHash[chainB];
                }
                else
                {
                    chainBSa = -1;
                }

                if (chainASa > -1 && chainBSa > -1)
                {
                    asa = interfaceAsa.CalculateInterfaceBuriedSurfaceArea(chainASa, chainBSa, complexAsa);
                    entryInterfaces[i].surfaceArea = asa;
                }
                else
                {
                    entryInterfaces[i].surfaceArea = -1;
                }
            }
        } 
     
        /// <summary>
        /// 
        /// </summary>
        /// <param name="symString"></param>
        /// <returns></returns>
        private string GetChainFromSymString (string symString)
        {
            string[] fields = symString.Split('_');
            return fields[0];
        }
        #endregion

        #region update residue numbering to author residue numbering
        /// <summary>
        /// 
        /// </summary>
        /// <param name="entryInterfaces"></param>
        private void UpdateEntryInterfaceResidueNumberingToAuth (InterfaceChains[] entryInterfaces)
        {
            for (int i = 0; i < entryInterfaces.Length; i ++)
            {
                for (int j = 0; j < entryInterfaces[i].chain1.Length; j ++)
                {
                    entryInterfaces[i].chain1[j].seqId = entryInterfaces[i].chain1[j].authSeqId;
                    entryInterfaces[i].chain1[j].residue = entryInterfaces[i].chain1[j].authResidue;
                }

                for (int j = 0; j < entryInterfaces[i].chain2.Length; j++)
                {
                    entryInterfaces[i].chain2[j].seqId = entryInterfaces[i].chain2[j].authSeqId;
                    entryInterfaces[i].chain2[j].residue = entryInterfaces[i].chain2[j].authResidue;
                }
                entryInterfaces[i].GetInterfaceResidueDist();
            }
        }
        #endregion
    }
}
