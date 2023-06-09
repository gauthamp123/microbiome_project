import xml.etree.ElementTree as ET
import xml

# GCF_009648975.1/xml/WP_064323585.1.xml
xmlfile = "GCF_009648975.1/xml/WP_064323585.1.xml"

tree = ET.parse(xmlfile)
root = tree.getroot()

print(root)

for item in root.findall('./BlastOutput_iterations/Iteration/Iteration_hits/Hit'):
    print(item)
    print(item.find('Hit_def'))