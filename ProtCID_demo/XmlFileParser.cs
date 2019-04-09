using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Net;
using System.IO;
using System.Data;
using System.Xml;
using AuxFuncLib;

namespace ProtCID_demo
{
    public class XmlFileParser
    {       
        private WebClient webClient = new WebClient();
        private const string PdbHttpDownloadAddress = "https://files.rcsb.org/download/";
        private string localXmlDir = "";

        public XmlFileParser (string xmlDir)
        {
            localXmlDir = xmlDir;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="xmlFile"></param>
        /// <returns></returns>
        public DataTable[] RetrieveEntryInfoFromXmlFile(string xmlFile)
        {           
            FileInfo fileInfo = new FileInfo(xmlFile);
            string pdbId = fileInfo.Name.Substring(0, 4);
            if (!File.Exists(xmlFile))
            {
                if (!File.Exists(xmlFile + ".gz"))
                {
                    DownloadXmlFile(fileInfo.Name);
                }
                else
                {
                    ParseHelper.UnZipFile(xmlFile + ".gz");
                }
            }

            if (xmlFile.IndexOf(".xml.gz") > -1)
            {
                ParseHelper.UnZipFile(xmlFile);
                xmlFile = xmlFile.Replace(".gz", "");
            }

            XmlDocument xmlDoc = new XmlDocument();
            xmlDoc.Load(xmlFile);
            // Create an XmlNamespaceManager for resolving namespaces.
            XmlNamespaceManager nsManager = new XmlNamespaceManager(xmlDoc.NameTable);
            string xmlNameSpace = xmlDoc.DocumentElement.Attributes["xmlns:PDBx"].InnerText;
            //		nsManager.AddNamespace("PDBx", "http://deposit.pdb.org/pdbML/pdbx.xsd");
            nsManager.AddNamespace("PDBx", xmlNameSpace);

            DataTable entityChainInfoTable = InitializeEntityChainInfoTable();            
            XmlNode polySchemeNode = xmlDoc.DocumentElement.SelectSingleNode("descendant::PDBx:pdbx_poly_seq_schemeCategory", nsManager);
            ParseEntityChainInfoPolyScheme(pdbId, polySchemeNode, entityChainInfoTable);

            XmlNode nonPolySchemeNode = xmlDoc.DocumentElement.SelectSingleNode("descendant::PDBx:pdbx_nonpoly_schemeCategory", nsManager);
            ParseNonPolyScheme(pdbId, nonPolySchemeNode, entityChainInfoTable);

            // get polymer type for each entity number
            // Polytype can be polypeptide, polydeoxyribonucleotide, polyribonucleotide or polysaccharide
            XmlNode entityPolyCategoryNode = xmlDoc.DocumentElement.SelectSingleNode("descendant::PDBx:entity_polyCategory", nsManager);
            ParseEntityPolyCategory(entityPolyCategoryNode, entityChainInfoTable);

            XmlNode nameCategoryNode = xmlDoc.DocumentElement.SelectSingleNode("descendant::PDBx:entity_name_comCategory", nsManager);
            ParseNameCategory(nameCategoryNode, entityChainInfoTable);

            DataTable entityDbInfoTable = InitializeEntityUnpTable();
            ParseDbRefCategory(pdbId, xmlDoc, nsManager, entityDbInfoTable);

            DataTable entityDbSeqTable = InitializeEntityDbSeqTable(); 
            ParseDbRefSeqCategory(pdbId, xmlDoc, nsManager, entityDbSeqTable);

            DataTable entryInfoTable = InitializeEntryInfoTable ();
            ParseEntryInfo(pdbId, xmlDoc, nsManager, entryInfoTable);

            DataTable[] entryInfoTables = new DataTable[4];
            entryInfoTables[0] = entityChainInfoTable;
            entryInfoTables[1] = entityDbInfoTable;
            entryInfoTables[2] = entryInfoTable;
            entryInfoTables[3] = entityDbSeqTable;
            return entryInfoTables;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        private DataTable InitializeEntityChainInfoTable ()
        {
            DataTable entryInfoTable = new DataTable("EntityChainInfo");
            string[] chainInfoColumns = { "PdbID", "AsymID", "EntityID", "AuthorChain", "Sequence", "CoordSequence",
                                            "PdbSeqNumbers", "AuthSeqNumbers", "DbSeqNumbers", "PolymerType", "PolymerStatus", "Name" };
            foreach (string col in chainInfoColumns)
            {
                entryInfoTable.Columns.Add(new DataColumn(col));
            }
            return entryInfoTable;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        private DataTable InitializeEntityUnpTable ()
        {
            DataTable entityUnpTable = new DataTable("EntityUnpDbRef");
            string[] unpInfoCols = {"PdbID", "RefID", "EntityID", "DbCode", "DBName", "DbAccession"};
            foreach (string unpCol in unpInfoCols)
            {
                entityUnpTable.Columns.Add(new DataColumn(unpCol));
            }
            return entityUnpTable;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        private DataTable InitializeEntityDbSeqTable()
        {
            DataTable entityDbSeqTable = new DataTable("EntityDbSeqRef");
            string[] dbRefInfoColumns = {"PdbID", "AlignID", "RefID", "DbAlignBeg", "DbAlignEnd", "AuthorChain", "AuthorAlignBeg", 
                                        "AuthorAlignEnd", "SeqAlignBeg", "SeqAlignEnd"};
            foreach (string col in dbRefInfoColumns)
            {
                entityDbSeqTable.Columns.Add(new DataColumn (col));
            }
            return entityDbSeqTable;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        private DataTable InitializeEntryInfoTable ()
        {
            DataTable entryInfoTable = new DataTable("EntryInfo");
            string[] entryInfoCols = {"PdbID", "Resolution", "Method", "Length_a", "Length_b", "Length_c", 
                                     "Angle_alpha", "Angle_beta", "Angle_gamma", "Z_PDB", "SpaceGroup", "Title", "Descript", 
                                     "CrystDetails", "CrystMethod", "Tempareture", "PH", "DepositFileDate", "ReleaseFileDate", "NumOfLigandAtoms", 
                                     "Rfactor_Rwork", "Rfactor_Rfree"};
            foreach (string colName in entryInfoCols)
            {
                entryInfoTable.Columns.Add(new DataColumn(colName));
            }
            return entryInfoTable;
        }

        /// <summary>
        /// get author/PDB asym correspondence for polymers
        /// </summary>
        /// <param name="polySchemeNode"></param>
        /// <param name="asymIdInfoHash"></param>
        private void ParseEntityChainInfoPolyScheme(string pdbId, XmlNode polySchemeNode, DataTable entityChainInfoTable)
        {
            if (polySchemeNode != null)
            {
                XmlNodeList asymPolyNodeList = polySchemeNode.ChildNodes;

                string preAsymId = "";
                string asymId = "";
                string entityId = "";
                string strandId = "";
                string residue = "";
                string resSequence = "";
                string resCoordinate = "";
                string authSeqNumbers = "";
                string pdbSeqNumbers = "";
                string ndbSeqNumbers = "";
                string insCode = "";

                if (asymPolyNodeList.Count > 0)
                {
                    preAsymId = asymPolyNodeList[0].Attributes["asym_id"].InnerText.ToString();
                    entityId = asymPolyNodeList[0].Attributes["entity_id"].InnerText.ToString();
                }

                foreach (XmlNode asymPolyNode in asymPolyNodeList)
                {
                    asymId = asymPolyNode.Attributes["asym_id"].InnerText.ToString();
                    if (preAsymId != asymId) // start new asymId
                    {
                        DataRow dataRow = entityChainInfoTable.NewRow();
                        dataRow["PdbID"] = pdbId;
                        dataRow["AsymID"] = preAsymId;
                        dataRow["EntityID"] = entityId;
                        dataRow["AuthorChain"] = strandId;
                        dataRow["Sequence"] = resSequence;
                        dataRow["CoordSequence"] = resCoordinate;
                        dataRow["AuthSeqNumbers"] = authSeqNumbers;
                        dataRow["PdbSeqNumbers"] = pdbSeqNumbers;
                        dataRow["DbSeqNumbers"] = ndbSeqNumbers;
                        dataRow["PolymerType"] = "poly";
                        dataRow["PolymerStatus"] = "polymer";
                        entityChainInfoTable.Rows.Add(dataRow);
                        // reset
                        preAsymId = asymId;
                        resSequence = "";
                        resCoordinate = "";
                        authSeqNumbers = "";
                        ndbSeqNumbers = "";
                        pdbSeqNumbers = "";
                        entityId = asymPolyNode.Attributes["entity_id"].InnerText.ToString();
                    }

                    residue = asymPolyNode.Attributes["mon_id"].InnerText.ToString();
                    if (residue.Length == 3)
                    {
                        string oneLetterResidue = ParseHelper.threeToOne(residue);
                        resSequence += oneLetterResidue;
                    }
                    else
                    {
                        resSequence += residue;
                    }

                    strandId = "_"; // defualt value is underscore

                    XmlNodeList polySeqNodeList = asymPolyNode.ChildNodes;
                    residue = "-";
                    insCode = "";
                    foreach (XmlNode polySeqNode in polySeqNodeList)
                    {
                        if (polySeqNode.Name.ToLower().IndexOf("strand_id") > -1)
                        {
                            strandId = polySeqNode.InnerText.ToString();
                            if (strandId.Length == 0)
                            {
                                strandId = "_";
                            }
                        }
                        if (polySeqNode.Name.ToLower().IndexOf("pdb_mon_id") > -1)
                        {
                            if (polySeqNode.InnerText != "")
                            {
                                residue = polySeqNode.InnerText;
                            }
                            if (residue.Length == 3)
                            {
                                residue = ParseHelper.threeToOne(residue);
                            }
                        }
                        if (polySeqNode.Name.ToLower().IndexOf("ndb_seq_num") > -1)
                        {
                            ndbSeqNumbers += polySeqNode.InnerText;
                        }
                        if (polySeqNode.Name.ToLower().IndexOf("pdb_seq_num") > -1)
                        {
                            pdbSeqNumbers += polySeqNode.InnerText;
                        }
                        if (polySeqNode.Name.ToLower().IndexOf("auth_seq_num") > -1)
                        {
                            authSeqNumbers += polySeqNode.InnerText;
                        }
                        if (polySeqNode.Name.ToLower().IndexOf("pdb_ins_code") > -1)
                        {
                            insCode = polySeqNode.InnerText;
                        }
                    }
                    ndbSeqNumbers += ",";
                    if (insCode != "")
                    {
                        pdbSeqNumbers += (insCode + ",");
                        authSeqNumbers += (insCode + ",");
                    }
                    else
                    {
                        pdbSeqNumbers += ",";
                        authSeqNumbers += ",";
                    }
                    resCoordinate += residue;
                }
                DataRow lastDataRow = entityChainInfoTable.NewRow();
                lastDataRow["PdbID"] = pdbId;
                lastDataRow["AsymID"] = preAsymId;
                lastDataRow["EntityID"] = entityId;
                lastDataRow["AuthorChain"] = strandId;
                lastDataRow["PolymerType"] = "poly";
                lastDataRow["PolymerStatus"] = "polymer";
                lastDataRow["Sequence"] = resSequence;
                lastDataRow["CoordSequence"] = resCoordinate;
                lastDataRow["AuthSeqNumbers"] = authSeqNumbers;
                lastDataRow["PdbSeqNumbers"] = pdbSeqNumbers;
                lastDataRow["DbSeqNumbers"] = ndbSeqNumbers;
                entityChainInfoTable.Rows.Add(lastDataRow);
            }
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="nonPolySchemeNode"></param>
        /// <param name="asymIdInfoHash"></param>
        /// <param name="strandIdAsymIdHash"></param>
        private void ParseNonPolyScheme(string pdbId, XmlNode nonPolySchemeNode, DataTable entityChainInfoTable)
        {
            if (nonPolySchemeNode != null)
            {
                string preAsymId = "";
                string asymId = "";
                string entityId = "";
                string strandId = "";
                string pdbMonId = "";
                string resSequence = "";
                string resCoordinate = "";
                string authSeqNumbers = "";
                string ndbSeqNumbers = "";
                string pdbSeqNumbers = "";

                XmlNodeList asymNonPolyNodeList = nonPolySchemeNode.ChildNodes;
                if (asymNonPolyNodeList.Count > 0)
                {
                    preAsymId = asymNonPolyNodeList[0].Attributes["asym_id"].InnerText.ToString();
                }

                foreach (XmlNode asymNonPolyNode in asymNonPolyNodeList)
                {
                    asymId = asymNonPolyNode.Attributes["asym_id"].InnerText.ToString();

                    if (preAsymId != asymId) // start new asymId
                    {
                        DataRow dataRow = entityChainInfoTable.NewRow();
                        dataRow["PdbID"] = pdbId;
                        dataRow["AsymID"] = preAsymId;
                        dataRow["EntityID"] = entityId;
                        dataRow["AuthorChain"] = strandId;
                        dataRow["Sequence"] = resSequence;
                        dataRow["CoordSequence"] = resCoordinate;
                        dataRow["AuthSeqNumbers"] = authSeqNumbers;
                        dataRow["PdbSeqNumbers"] = pdbSeqNumbers;
                        dataRow["DbSeqNumbers"] = ndbSeqNumbers;
                        dataRow["PolymerType"] = "nonpoly";
                        dataRow["PolymerStatus"] = "non-polymer";
                        entityChainInfoTable.Rows.Add(dataRow);
                        // reset
                        preAsymId = asymId;
                        resSequence = "";
                        resCoordinate = "";
                        authSeqNumbers = "";
                        ndbSeqNumbers = "";
                        pdbSeqNumbers = "";
                    }

                    strandId = "_"; // defualt value is underscore

                    ndbSeqNumbers += asymNonPolyNode.Attributes["ndb_seq_num"].InnerText.ToString();

                    XmlNodeList nonpolyNodeList = asymNonPolyNode.ChildNodes;
                    pdbMonId = "-";
                    foreach (XmlNode nonpolyNode in nonpolyNodeList)
                    {
                        if (nonpolyNode.Name.ToLower().IndexOf("entity_id") > -1)
                        {
                            entityId = nonpolyNode.InnerText.ToString();
                        }
                        if (nonpolyNode.Name.ToLower().IndexOf("strand_id") > -1)
                        {
                            strandId = nonpolyNode.InnerText.ToString();
                            if (strandId.Length == 0)
                                strandId = "_";
                        }
                        if (nonpolyNode.Name.ToLower().IndexOf("pdbx:mon_id") > -1)
                        {
                            resSequence += nonpolyNode.InnerText;
                        }
                        if (nonpolyNode.Name.ToLower().IndexOf("pdb_mon_id") > -1)
                        {
                            pdbMonId = nonpolyNode.InnerText;
                        }
                        if (nonpolyNode.Name.ToLower().IndexOf("pdb_seq_num") > -1)
                        {
                            pdbSeqNumbers += nonpolyNode.InnerText;
                        }
                        if (nonpolyNode.Name.ToLower().IndexOf("auth_seq_num") > -1)
                        {
                            authSeqNumbers += nonpolyNode.InnerText;
                        }
                    }
                    authSeqNumbers += ",";
                    ndbSeqNumbers += ",";
                    pdbSeqNumbers += ",";
                    resCoordinate += pdbMonId;
                }
                // add the last asymid
                DataRow lastDataRow = entityChainInfoTable.NewRow();
                lastDataRow["PdbID"] = pdbId;
                lastDataRow["AsymID"] = preAsymId;
                lastDataRow["EntityID"] = entityId;
                lastDataRow["AuthorChain"] = strandId;
                lastDataRow["Sequence"] = resSequence;
                lastDataRow["CoordSequence"] = resCoordinate;
                lastDataRow["AuthSeqNumbers"] = authSeqNumbers;
                lastDataRow["PdbSeqNumbers"] = pdbSeqNumbers;
                lastDataRow["DbSeqNumbers"] = ndbSeqNumbers;
                lastDataRow["PolymerType"] = "nonpoly";
                lastDataRow["PolymerStatus"] = "non-polymer";
                entityChainInfoTable.Rows.Add(lastDataRow);
            }
        }

        /// <summary>
        /// get polymer type for each entity number
        /// </summary>
        /// <param name="entityPolyCategoryNode"></param>
        /// <param name="entityIdPolyTypeHash"></param>
        private void ParseEntityPolyCategory(XmlNode entityPolyCategoryNode, DataTable entityChainInfoTable)
        {
            if (entityPolyCategoryNode != null)
            {
                XmlNodeList entityPolyNodeList = entityPolyCategoryNode.ChildNodes;
                string entityPolyType = "";
                foreach (XmlNode entityPolyNode in entityPolyNodeList)
                {
                    entityPolyType = "";

                    string entityId = entityPolyNode.Attributes[0].InnerText;
                    foreach (XmlNode polyTypeNode in entityPolyNode.ChildNodes)
                    {
                        if (polyTypeNode.Name.ToLower() == "pdbx:type")
                        {
                            entityPolyType = polyTypeNode.InnerText;
                        }                       
                    }
                    if (entityPolyType.Length == 0)
                    {
                        entityPolyType = "-";
                    }
                    if (entityPolyType.IndexOf ("polypeptide") > -1)
                    {
                        entityPolyType = "polypeptide";
                    }
                    DataRow[] entityRows = entityChainInfoTable.Select(string.Format ("EntityID = '{0}'", entityId));
                    foreach (DataRow entityRow in entityRows)
                    {
                        entityRow["PolymerType"] = entityPolyType;
                    }
                }
            }
        }

        /// <summary>
        /// get name source for entity id
        /// </summary>
        /// <param name="nameCategoryNode"></param>
        /// <param name="entityIdRcsbNameHash"></param>
        /// <param name="entityIdSwsHash"></param>
        private void ParseNameCategory(XmlNode nameCategoryNode, DataTable entityChainInfoTable)
        {
            if (nameCategoryNode != null)
            {
                XmlNodeList entityNameNodeList = nameCategoryNode.ChildNodes;
                string entityName = "";
                foreach (XmlNode entityNameNode in entityNameNodeList)
                {
                    string entityId = entityNameNode.Attributes["entity_id"].InnerText;
                    entityName = "";
                    if (entityNameNode["PDBx:name"] != null)
                    {
                        entityName = entityNameNode["PDBx:name"].InnerText.ToString().ToUpper();
                    }
                    DataRow[] entityRows = entityChainInfoTable.Select(string.Format ("EntityID = '{0}'", entityId));
                    foreach (DataRow entityRow in entityRows)
                    {
                        entityRow["Name"] = entityName;
                    }
                }
            }
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="pdbId"></param>
        /// <param name="xmlDoc"></param>
        /// <param name="nsManager"></param>
        private void ParseDbRefCategory(string pdbId, XmlDocument xmlDoc, XmlNamespaceManager nsManager, DataTable entityUnpCodeTable)
        {
            XmlNode dbRefCatNode = xmlDoc.DocumentElement.SelectSingleNode("descendant::PDBx:struct_refCategory", nsManager);
            if (dbRefCatNode == null)
            {
                return;
            }
            string dbCode = "";
            string dbName = "";
            string dbAccession = "";
            int entityId = 0;
            int refId = 0;
            foreach (XmlNode dbRefNode in dbRefCatNode.ChildNodes)
            {
                entityId = Convert.ToInt32 (dbRefNode["PDBx:entity_id"].InnerText);
                refId = Convert.ToInt32 (dbRefNode.Attributes["id"].InnerText);
                if (dbRefNode["PDBx:db_code"] == null || dbRefNode["PDBx:db_code"].InnerText == "")
                {
                    dbCode = "-";
                }
                else
                {
                    dbCode = dbRefNode["PDBx:db_code"].InnerText;
                }
                if (dbRefNode["PDBx:db_name"] == null || dbRefNode["PDBx:db_name"].InnerText == "")
                {
                   dbName = "-";
                }
                else
                {
                    dbName = dbRefNode["PDBx:db_name"].InnerText;
                }
                if (dbRefNode["PDBx:pdbx_db_accession"] == null || dbRefNode["PDBx:pdbx_db_accession"].InnerText == "")
                {
                    dbAccession = "-";                    
                }
                else
                {
                    dbAccession = dbRefNode["PDBx:pdbx_db_accession"].InnerText;
                }
                if (dbName == "UNP")
                {
                    DataRow unpRow = entityUnpCodeTable.NewRow();
                    unpRow["PdbID"] = pdbId;
                    unpRow["EntityID"] = entityId;
                    unpRow["RefID"] = refId;
                    unpRow["DbName"] = dbName;
                    unpRow["DBCode"] = dbCode;
                    unpRow["DBAccession"] = dbAccession;
                    entityUnpCodeTable.Rows.Add(unpRow);
                }
            }
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="pdbId"></param>
        /// <param name="xmlDoc"></param>
        /// <param name="nsManager"></param>
        private void ParseDbRefSeqCategory(string pdbId, XmlDocument xmlDoc, XmlNamespaceManager nsManager, DataTable entityDbSeqTable)
        {
            XmlNode dbRefSeqCatNode = xmlDoc.DocumentElement.SelectSingleNode("descendant::PDBx:struct_ref_seqCategory", nsManager);
            if (dbRefSeqCatNode == null)
            {
                Console.WriteLine("No PDBx:struct_ref_seqCategory exist in " + pdbId);
                return;
            }
            foreach (XmlNode dbRefSeqNode in dbRefSeqCatNode.ChildNodes)
            {
                DataRow dbRefSeqRow = entityDbSeqTable.NewRow();
                dbRefSeqRow["PdbID"] = pdbId;
                dbRefSeqRow["AlignID"] = dbRefSeqNode.Attributes["align_id"].InnerText;
                dbRefSeqRow["RefID"] = dbRefSeqNode["PDBx:ref_id"].InnerText;
                dbRefSeqRow["DbAlignBeg"] = dbRefSeqNode["PDBx:db_align_beg"].InnerText;
                dbRefSeqRow["DbAlignEnd"] = dbRefSeqNode["PDBx:db_align_end"].InnerText;
                if (dbRefSeqNode["PDBx:pdbx_strand_id"] != null)
                {
                    dbRefSeqRow["AuthorChain"] = dbRefSeqNode["PDBx:pdbx_strand_id"].InnerText;
                }
                else
                {
                    dbRefSeqRow["AuthorChain"] = "-";
                }
                if (dbRefSeqNode["PDBx:pdbx_auth_seq_align_beg"] != null)
                {
                    dbRefSeqRow["AuthorAlignBeg"] = dbRefSeqNode["PDBx:pdbx_auth_seq_align_beg"].InnerText;
                }
                else
                {
                    dbRefSeqRow["AuthorAlignBeg"] = -1;
                }
                if (dbRefSeqNode["PDBx:pdbx_auth_seq_align_end"] != null)
                {
                    dbRefSeqRow["AuthorAlignEnd"] = dbRefSeqNode["PDBx:pdbx_auth_seq_align_end"].InnerText;
                }
                else
                {
                    dbRefSeqRow["AuthorAlignEnd"] = -1;
                }
                dbRefSeqRow["SeqAlignBeg"] = dbRefSeqNode["PDBx:seq_align_beg"].InnerText;
                dbRefSeqRow["SeqAlignEnd"] = dbRefSeqNode["PDBx:seq_align_end"].InnerText;
                entityDbSeqTable.Rows.Add(dbRefSeqRow);
            }
        }

        #region parse entry info from xml file
        /// <summary>
        /// 
        /// </summary>
        /// <param name="pdbId"></param>
        /// <param name="authChain"></param>
        /// <returns></returns>
        private int GetAuthChainEntityIdFromXmlFile(string pdbId, string authChain)
        {
            string xmlFile = Path.Combine(localXmlDir, pdbId + ".xml.gz");
            if (!File.Exists(xmlFile))
            {
                DownloadXmlFile(pdbId);
            }
            int entityId = 0;
            return entityId;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="xmlDoc"></param>
        /// <param name="entryInfo"></param>
        /// <param name="nsManager"></param>
        private void ParseEntryInfo(string pdbId, XmlDocument xmlDoc, XmlNamespaceManager nsManager, DataTable entryInfoTable)
        {
            DataRow entryInfoRow = entryInfoTable.NewRow();
            entryInfoRow["PdbID"] = pdbId;
            // get number of atoms ligand
            XmlNode numberLigandAtomsNode = xmlDoc.DocumentElement.SelectSingleNode("descendant::PDBx:refine_histCategory/PDBx:refine_hist/PDBx:pdbx_number_atoms_ligand", nsManager);
            if (numberLigandAtomsNode == null)
            {
                entryInfoRow["NumOfLigandAtoms"] = 0;
            }
            else
            {
                string numOfLigandString = numberLigandAtomsNode.InnerText.ToString();
                if (numOfLigandString == "")
                {
                    entryInfoRow["NumOfLigandAtoms"] = 0;
                }
                else
                {
                    entryInfoRow["NumOfLigandAtoms"] = System.Convert.ToInt32(numOfLigandString);
                }
            }

            // get resolution from refineCategory
            XmlNode resolutionNode = xmlDoc.DocumentElement.SelectSingleNode("descendant::PDBx:refineCategory/PDBx:refine/PDBx:ls_d_res_high", nsManager);
            if (resolutionNode == null) // if EM structures, added on Jan. 3, 2019
            {
                resolutionNode = xmlDoc.DocumentElement.SelectSingleNode("descendant::PDBx:em_3d_reconstructionCategory/PDBx:em_3d_reconstruction/PDBx:resolution", nsManager);
            }

            if (resolutionNode != null)
            {
                string resStr = resolutionNode.InnerText.ToString();
                if (resStr != "")
                {
                    entryInfoRow["Resolution"] = System.Convert.ToDouble(resStr);
                }
            }

            // get R factor R work from refineCategory
            XmlNode rfactorworkNode = xmlDoc.DocumentElement.SelectSingleNode("descendant::PDBx:refineCategory/PDBx:refine/PDBx:ls_R_factor_R_work", nsManager);
            if (rfactorworkNode != null)
            {
                string resStr = rfactorworkNode.InnerText.ToString();
                if (resStr != "")
                {
                    entryInfoRow["Rfactor_Rwork"] = System.Convert.ToDouble(resStr);
                }
            }

            // get R Factor R Free from refineCategory
            XmlNode rfactorfreeNode = xmlDoc.DocumentElement.SelectSingleNode("descendant::PDBx:refineCategory/PDBx:refine/PDBx:ls_R_factor_R_free", nsManager);
            if (rfactorfreeNode != null)
            {
                string resStr = rfactorfreeNode.InnerText.ToString();
                if (resStr != "")
                {
                    entryInfoRow["Rfactor_Rfree"] = System.Convert.ToDouble(resStr);
                }
            }

            // get method
            XmlNode methodNode = xmlDoc.DocumentElement.SelectSingleNode("descendant::PDBx:exptlCategory/PDBx:exptl", nsManager);
            if (methodNode == null)
            {
                entryInfoRow["Method"] = "-";
            }
            else
            {
                if (methodNode["PDBx:method"] != null)
                {
                    entryInfoRow["Method"] = methodNode["PDBx:method"].InnerText;
                }
                else
                {
                    entryInfoRow["Method"] = methodNode.Attributes["method"].Value;
                }
            }

            // get cell information: Length, Angle, and Z_PDB
            string cellInfoNodePath = "descendant::PDBx:cellCategory/PDBx:cell";
            XmlNode cellInfoNode = xmlDoc.DocumentElement.SelectSingleNode(cellInfoNodePath, nsManager);
            if (cellInfoNode != null)
            {
                XmlNodeList cellInfoNodeList = cellInfoNode.ChildNodes;
                double cellValue = -1.0;
                foreach (XmlNode cellNode in cellInfoNodeList)
                {
                    string nodeName = cellNode.Name.ToLower();
                    if (cellNode.InnerText.ToString() == "")
                    {
                        cellValue = -1.0;
                    }
                    else
                    {
                        cellValue = System.Convert.ToDouble(cellNode.InnerText.ToString());
                    }
                    switch (nodeName)
                    {
                        case "pdbx:length_a":
                            entryInfoRow["Length_a"] = cellValue;
                            break;
                        case "pdbx:length_b":
                            entryInfoRow["Length_b"] = cellValue;
                            break;
                        case "pdbx:length_c":
                            entryInfoRow["Length_c"]= cellValue;
                            break;
                        case "pdbx:angle_alpha":
                            entryInfoRow["Angle_alpha"] = cellValue;
                            break;
                        case "pdbx:angle_beta":
                            entryInfoRow["Angle_beta"] = cellValue;
                            break;
                        case "pdbx:angle_gamma":
                            entryInfoRow["Angle_gamma"] = cellValue;
                            break;
                        case "pdbx:z_pdb":
                            entryInfoRow["Z_PDB"] = (int)cellValue;
                            break;
                    }
                }
            }

            XmlNode spaceGroupNode = xmlDoc.DocumentElement.SelectSingleNode("descendant::PDBx:symmetryCategory/PDBx:symmetry/PDBx:space_group_name_H-M", nsManager);
            if (spaceGroupNode == null)
            {
                entryInfoRow["SpaceGroup"] = "-";
            }
            else
            {
                entryInfoRow["SpaceGroup"] = spaceGroupNode.InnerText.ToString();
            }

            // get title and description for the entry
            XmlNode structNode = xmlDoc.DocumentElement.SelectSingleNode("descendant::PDBx:structCategory/PDBx:struct", nsManager);
            if (structNode == null)
            {
                entryInfoRow["Title"] = "-";
                entryInfoRow["Descript"] = "-";
            }
            else
            {
                XmlNodeList childNodes = structNode.ChildNodes;
                foreach (XmlNode cellNode in childNodes)
                {
                    string nodeName = cellNode.Name.ToLower();
                    switch (nodeName)
                    {
                        case "pdbx:title":
                            entryInfoRow["Title"] = cellNode.InnerText.ToString().ToUpper();
                            break;

                        case "pdbx:pdbx_descriptor":
                            entryInfoRow["Descript"] = cellNode.InnerText.ToString().ToUpper();
                            break;

                        default:
                            break;
                    }
                }
            }
            // get crystalization conditions for the crystal structure
            XmlNode crystConditionNode = xmlDoc.DocumentElement.SelectSingleNode("descendant::PDBx:exptl_crystal_growCategory/PDBx:exptl_crystal_grow", nsManager);
            entryInfoRow["PH"] = -1;
            entryInfoRow["Tempareture"]= -1;
            entryInfoRow["CrystDetails"] = "";
            entryInfoRow["CrystMethod"] = "";

            if (crystConditionNode != null)
            {
                string tempString = "";
                string phString = "";
                XmlNodeList childNodes = crystConditionNode.ChildNodes;
                foreach (XmlNode crystNode in childNodes)
                {
                    string nodeName = crystNode.Name.ToLower();
                    switch (nodeName)
                    {
                        case "pdbx:temp":
                            tempString = crystNode.InnerText.ToString().TrimEnd('.');
                            if (tempString != "")
                            {
                                entryInfoRow["Tempareture"] = Convert.ToDouble(tempString);
                            }
                            else
                            {
                                entryInfoRow["Tempareture"] = -1.0;
                            }
                            break;

                        case "pdbx:ph":
                            phString = crystNode.InnerText.ToString();
                            if (phString != "")
                            {
                                entryInfoRow["PH"] = Convert.ToDouble(phString);
                            }
                            else
                            {
                                entryInfoRow["PH"] = -1.0;
                            }
                            break;

                        case "pdbx:pdbx_details":
                            entryInfoRow["CrystDetails"] = crystNode.InnerText.ToString();
                            break;

                        case "pdbx:method":
                            entryInfoRow["CrystMethod"] = crystNode.InnerText.ToString();
                            break;

                        default:
                            break;
                    }
                }
            }
            entryInfoTable.Rows.Add(entryInfoRow);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="pdbId"></param>
        private void DownloadXmlFile(string xmlFileName)
        {
            string httpFile = PdbHttpDownloadAddress + xmlFileName;
            string localFile = Path.Combine(localXmlDir, xmlFileName);
            webClient.DownloadFile(httpFile, localFile);
        }
        #endregion

    }
}
