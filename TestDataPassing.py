import unittest
import DataTreatment 
from DataTreatment import getData, openJson, DataType

class CopyData(unittest.TestCase):
    '''Class to test how data is copied between functions.'''
    expected_simple_return = [
        {
        "wild-type": False, 
        "organism": "Mus musculus", 
        "turnover": "0.02", 
    }, 
    {
        "wild-type": True, 
        "organism": "Caenorhabditis elegans", 
        "turnover": "0.052", 
    }, 
    {
        "wild-type": False, 
        "organism": "Rattus norvegicus", 
        "turnover": "0.243", 
    }, 
    {
        "wild-type": False, 
        "organism": "Solanum lycopersicum", 
        "turnover": "0.38", 
    }, 
    {
        "wild-type": False, 
        "organism": "Entamoeba histolytica", 
        "turnover": "0.5", 
    }, 
    {
        "wild-type": False, 
        "organism": "Entamoeba histolytica", 
        "turnover": "1.25", 
    }, 
    {
        "wild-type": False, 
        "organism": "Sulfolobus solfataricus", 
        "turnover": "1.3", 
    }, 
    {
        "wild-type": False, 
        "organism": "Escherichia coli", 
        "turnover": "2.38", 
    }, 
    {
        "wild-type": False, 
        "organism": "Mus musculus", 
        "turnover": "3.5", 
    }, 
    {
        "wild-type": False, 
        "organism": "Caenorhabditis elegans", 
        "turnover": "4.03", 
    }, 
    {
        "wild-type": False, 
        "organism": "Caenorhabditis elegans", 
        "turnover": "5.3", 
    }, 
    {
        "wild-type": False, 
        "organism": "Homo sapiens", 
        "turnover": "5.58", 
    }, 
    {
        "wild-type": False, 
        "organism": "Homo sapiens", 
        "turnover": "8.1", 
    }, 
    {
        "wild-type": True, 
        "organism": "Caenorhabditis elegans", 
        "turnover": "10.17", 
    }, 
    {
        "wild-type": False, 
        "organism": "Escherichia coli", 
        "turnover": "13.2", 
    }, 
    {
        "wild-type": False, 
        "organism": "Anopheles gambiae", 
        "turnover": "14.3", 
    }, 
    {
        "wild-type": False, 
        "organism": "Anopheles gambiae", 
        "turnover": "15.4", 
    }, 
    {
        "wild-type": False, 
        "organism": "Anopheles gambiae", 
        "turnover": "15.7", 
    }, 
    {
        "wild-type": False, 
        "organism": "Mus musculus", 
        "turnover": "19.97", 
    }, 
    {
        "wild-type": True, 
        "organism": "Escherichia coli", 
        "turnover": "22", 
    }, 
    {
        "wild-type": True, 
        "organism": "Escherichia coli", 
        "turnover": "22.8", 
    }, 
    {
        "wild-type": True, 
        "organism": "Mus musculus", 
        "turnover": "25", 
    }, 
    {
        "wild-type": True, 
        "organism": "Homo sapiens", 
        "turnover": "25.78", 
    }, 
    {
        "wild-type": True, 
        "organism": "Homo sapiens", 
        "turnover": "27.4", 
    }, 
    {
        "wild-type": True, 
        "organism": "Mus musculus", 
        "turnover": "29.5", 
    }, 
    {
        "wild-type": True, 
        "organism": "Mus musculus", 
        "turnover": "37", 
    }, 
    {
        "wild-type": True, 
        "organism": "Mus musculus", 
        "turnover": "37.88", 
    }, 
    {
        "wild-type": True, 
        "organism": "Rattus norvegicus", 
        "turnover": "41.7", 
    }, 
    {
        "wild-type": True, 
        "organism": "Homo sapiens", 
        "turnover": "46.57", 
    }, 
    {
        "wild-type": True, 
        "organism": "Rattus norvegicus", 
        "turnover": "50", 
    }, 
    {
        "wild-type": True, 
        "organism": "Aeropyrum pernix", 
        "turnover": "63.2", 
    }, 
    {
        "wild-type": True, 
        "organism": "Bos taurus", 
        "turnover": "1030", 
    }, 
    {
        "wild-type": True, 
        "organism": "Bos taurus", 
        "turnover": "1200", 
    }, 
    {
        "wild-type": True, 
        "organism": "Bos taurus", 
        "turnover": "1300", 
    } ] 
    
    expected_turnover_return = [
        {
        "wild-type": False, 
        "organism": "Mus musculus", 
        "turnover": "0.02", 
    }, 
    {
        "wild-type": True, 
        "organism": "Caenorhabditis elegans", 
        "turnover": "0.052", 
    }, 
    {
        "wild-type": False, 
        "organism": "Rattus norvegicus", 
        "turnover": "0.243", 
    }, 
    {
        "wild-type": False, 
        "organism": "Solanum lycopersicum", 
        "turnover": "0.38", 
    }, 
    {
        "wild-type": False, 
        "organism": "Entamoeba histolytica", 
        "turnover": "0.5", 
    }, 
    {
        "wild-type": False, 
        "organism": "Entamoeba histolytica", 
        "turnover": "1.25", 
    }, 
    {
        "wild-type": False, 
        "organism": "Sulfolobus solfataricus", 
        "turnover": "1.3", 
    }]
    
    expected_specific_activity_return = [
        {
        "wild-type": False, 
        "organism": "Mus musculus", 
        "specific activity": "0.02", 
    }, 
    {
        "wild-type": True, 
        "organism": "Caenorhabditis elegans", 
        "specific activity": "0.052", 
    }, 
    {
        "wild-type": False, 
        "organism": "Rattus norvegicus", 
        "specific activity": "0.243", 
    }, 
    {
        "wild-type": False, 
        "organism": "Solanum lycopersicum", 
        "specific activity": "0.38", 
    }, 
    {
        "wild-type": False, 
        "organism": "Entamoeba histolytica", 
        "specific activity": "0.5", 
    }, 
    {
        "wild-type": False, 
        "organism": "Entamoeba histolytica", 
        "specific activity": "1.25", 
    }, 
    {
        "wild-type": False, 
        "organism": "Sulfolobus solfataricus", 
        "specific activity": "1.3", 
    }]
    
    expected_molecular_weight_return = [
        {
        "wild-type": False, 
        "organism": "Mus musculus", 
        "molecular weight": "0.02", 
    }, 
    {
        "wild-type": True, 
        "organism": "Caenorhabditis elegans", 
        "molecular weight": "0.052", 
    }, 
    {
        "wild-type": False, 
        "organism": "Rattus norvegicus", 
        "molecular weight": "0.243", 
    }, 
    {
        "wild-type": False, 
        "organism": "Solanum lycopersicum", 
        "molecular weight": "0.38", 
    }, 
    {
        "wild-type": False, 
        "organism": "Entamoeba histolytica", 
        "molecular weight": "0.5", 
    }, 
    {
        "wild-type": False, 
        "organism": "Entamoeba histolytica", 
        "molecular weight": "1.25", 
    }, 
    {
        "wild-type": False, 
        "organism": "Sulfolobus solfataricus", 
        "molecular weight": "1.3", 
    }]
        
    def test_getData(self):
        '''Ensure getData can extract and return correct data as a list of dict'''
        sample_reaction = openJson('Unit Tests/sample_brenda_reaction.json')
 
                
        
        returned_data = getData(sample_reaction, 'thioredoxin',  
                                DataType.turnover)
        for index, entry in enumerate(returned_data):
            self.assertDictEqual(returned_data[index], 
                                 CopyData.expected_simple_return[index])
     
    def test_get_turnover_data_from_mixed(self):
        '''Ensure getData can extract and return correct data as a list of dict
            when the data_types are mixed.'''
        sample_reaction = openJson('Unit Tests/sample_mixed_brenda_reaction.json')      
        returned_data = getData(sample_reaction, 'thioredoxin',  
                                DataType.turnover)
        for index, entry in enumerate(returned_data):
            self.assertDictEqual(returned_data[index], 
                                 CopyData.expected_turnover_return[index])
        
    def test_get_specific_activity_data_from_mixed(self):
        '''Ensure getData can extract specific activities data as a list of dict'''
        sample_reaction = openJson('Unit Tests/sample_mixed_brenda_reaction.json')      
        returned_data = getData(sample_reaction, 'thioredoxin',  
                                DataType.specific_activity)
        for index, entry in enumerate(returned_data):
            self.assertDictEqual(returned_data[index], 
                                 CopyData.expected_specific_activity_return[index])
            
    def test_get_molecular_weight_data_from_mixed(self):
        '''Ensure getData can extract specific activities data as a list of dict'''
        sample_reaction = openJson('Unit Tests/sample_mixed_brenda_reaction.json')      
        returned_data = getData(sample_reaction, 'thioredoxin',  
                                DataType.molecular_weight)
        for index, entry in enumerate(returned_data):
            self.assertDictEqual(returned_data[index], 
                                 CopyData.expected_molecular_weight_return[index])
        
if __name__ == '__main__':
    unittest.main()
    
