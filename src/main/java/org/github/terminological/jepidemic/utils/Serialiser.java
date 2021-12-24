package org.github.terminological.jepidemic.utils;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Paths;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import uk.co.terminological.rjava.RClass;
import uk.co.terminological.rjava.RMethod;
import uk.co.terminological.rjava.types.RDataframe;
import uk.co.terminological.rjava.types.RList;
import uk.co.terminological.rjava.types.RNamedList;
import uk.co.terminological.rjava.types.RObject;

/** Actually useful class which takes an R dataframe and a filename
 * string and serialises the dataframe so it can be used natively in java
 * for testing purposes.
 * @author terminological
 *
 */
@RClass
public class Serialiser {

	static Logger log = LoggerFactory.getLogger(Serialiser.class);
	
	//SNIPPET_1
	
	@RMethod
	public static void serialiseDataframe(RDataframe dataframe, String filename) throws IOException {
		FileOutputStream fos = new FileOutputStream(filename);
		dataframe.writeRDS(fos);
		log.info("dataframe written to: "+filename);
	}
	
	@RMethod
	public static RDataframe deserialiseDataframe(String filename) throws IOException {
		InputStream is = Files.newInputStream(Paths.get(filename));
		if(is==null) throw new IOException("Could not locate "+filename);
		return RObject.readRDS(RDataframe.class, is);
	}
	
	//SNIPPET_2
	
	@RMethod
	public static void serialiseList(RList dataframe, String filename) throws IOException {
		FileOutputStream fos = new FileOutputStream(filename);
		dataframe.writeRDS(fos);
		log.info("list written to: "+filename);
	}
	
	@RMethod
	public static RList deserialiseList(String filename) throws IOException {
		InputStream is = Files.newInputStream(Paths.get(filename));
		if(is==null) throw new IOException("Could not locate "+filename);
		return RObject.readRDS(RList.class, is);
	}
	
	//SNIPPET_3
	
	@RMethod
	public static void serialiseNamedList(RNamedList dataframe, String filename) throws IOException {
		FileOutputStream fos = new FileOutputStream(filename);
		dataframe.writeRDS(fos);
		log.info("named list written to: "+filename);
	}
	
	@RMethod
	public static RNamedList deserialiseNamedList(String filename) throws IOException {
		InputStream is = Files.newInputStream(Paths.get(filename));
		if(is==null) throw new IOException("Could not locate "+filename);
		return RObject.readRDS(RNamedList.class, is);
	}
	
	//SNIPPET_4
	
//	TODO: We don't support arbitary content unless it is a list, but we probably should for just this kind of situtation
//	@RMethod
//	public static void serialise(RObject anything, String filename) throws IOException {
//		FileOutputStream fos = new FileOutputStream(filename);
//		anything.writeRDS(fos);
//		log.info("named list written to: "+filename);
//	}
}
