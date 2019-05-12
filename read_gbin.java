//import gaia.cu3.avu.gsr.util.fileIO.FileIO;
import gaia.cu3.avu.gsr.util.fileIO.FileIO;
import gaia.cu3.avu.gsr.util.MultiUtility;
import gaia.cu3.dpct.avu.gsr.dm.GsrAstrometricSource;

public class read_gbin{
    public static void main( String[] args )
    {
    //ArrayList<GsrAstrometricSource> sources = loadAllFilesInDir(FileExtensionsEnum.GBIN);
	ArrayList<GsrAstrometricSource> sources = (ArrayList<GsrAstrometricSource>)loadGBinFile("GsrAstrometricSource_1597585_0000.gbin");
    	System.out.println(sources);
//    for(GsrAstrometricSource s : sources) {

//    	s.getDelta();

//    }
	}
}
