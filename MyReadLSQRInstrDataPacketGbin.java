//Utility class for read a GenCoeffDataPacket gBin
package gaia.cu3.avu.gsr.util.DataReader;

import gaia.cu1.tools.dal.ObjectFactory;
import gaia.cu1.tools.dal.file.GaiaRootGaiaTable;
import gaia.cu1.tools.exception.GaiaDataAccessException;
import gaia.cu1.tools.exception.GaiaException;
import gaia.cu1.tools.util.props.PropertyLoader;
import gaia.cu3.dpct.avu.bam.dm.BamBananaWrapper;
import gaia.cu3.dpct.avu.gsr.packet.dm.GenCoeffDataPacket;
import gaia.cu3.dpct.avu.gsr.packet.dm.InstrSegmentDataPacket;
import gaia.cu3.dpct.avu.gsr.packet.dm.LSQRInstrDataPacket;

import java.io.File;
import java.io.IOException;

/*public class MyReadBananaGbin <E extends GaiaRoot>extends GbinReaderAVU <E > {
	public MyReadBananaGbin(File file) throws GaiaDataAccessException {
		super(file);
		
	}*/
public class MyReadLSQRInstrDataPacketGbin  {
	File gbinIN;
	public MyReadLSQRInstrDataPacketGbin(File gbinIN) throws GaiaDataAccessException {
		
		this.gbinIN=gbinIN;	
	}
public LSQRInstrDataPacket [] readlSQRInstrDataPacketGbin() throws GaiaException, IOException {
	
		PropertyLoader.load();
//	     PropertyLoader.load( "conf/.properties", true );
		GaiaRootGaiaTable tableForLSQRInstrDataPacket = new GaiaRootGaiaTable();
		tableForLSQRInstrDataPacket.populateFromFile( gbinIN.getCanonicalPath() );

		ObjectFactory<LSQRInstrDataPacket> factoryForInst = new ObjectFactory<>( LSQRInstrDataPacket.class );
		LSQRInstrDataPacket[] packetArr = factoryForInst.getObjects( tableForLSQRInstrDataPacket );

		return packetArr;
		}

	
}
