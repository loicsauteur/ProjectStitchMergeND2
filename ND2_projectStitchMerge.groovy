/*
        Project, stitch and merge channels for nd2 files

Description:
    Allows you to project and stitch multipoint nd2 files.
    If you have several channels saved as individual files,
    this script will merge these channels. You will be prompted
    about the channels.

More details:
    Stitching of the projections is done in memory (no temporary files)
    Standard max projection for fluorescence and min projection for BF (only one channel)
    Background subtraction with 50px rolling ball, optional (on stacks)
    Channel information (if more than 1) is automatically identified for
        standard channels naming: DAPI, BF, 470, 555, 640
    No tile configuration file written
    nd2 files stage position for x is inverted

@author: Loïc Sauteur - loic.sauteur@unibas.ch - DBM, Basel

*/


// imports
import fiji.util.gui.GenericDialogPlus
import groovy.io.FileType
import ij.*
import ij.plugin.RGBStackMerge
import ij.plugin.ZProjector
import ij.plugin.filter.BackgroundSubtracter
import ij.process.LUT
import loci.formats.ImageReader
import loci.plugins.BF
import loci.plugins.in.ImporterOptions
import mpicbg.stitching.CollectionStitchingImgLib
import mpicbg.stitching.ImageCollectionElement
import mpicbg.stitching.StitchingParameters
import mpicbg.stitching.fusion.Fusion
import net.imglib2.img.display.imagej.ImageJFunctions
import mpicbg.models.TranslationModel2D

import java.awt.Color


// script parameters
#@ File(label="Choose a folder", style="directory") input
#@ File (label = "Output largeImages", style = "directory") output
#@ Integer (label = "Number of channels", value=4) nCh
#@ Boolean(label="Do a background correction?", value=false) bgCorrection

// global variables
suffix = ".nd2"

// script start
main()
println("Processing finished")


def main() {
    // create a file list (full path name of the files) of nd2 files in input
    // FIXME ! Does it also for sub directories
    def fileList = []
    def dirFile = new File(input.getAbsolutePath())
    dirFile.eachFileRecurse(FileType.FILES) { file ->
        if (file.name.endsWith(suffix)) {
            fileList << file
        }
    }
    fileList.sort()
    if (fileList.isEmpty()) {
        println("No .nd2 files found.")
        return
    }

    // identify the channels
    (chSeq, BFchNr) = getChannelIdentifier(fileList, nCh)

    // iterate over fileList    ********************************************************
    wellcounter = 1
    for (int i = 0; i < fileList.size() / nCh; i++) {
        println("Processing 'well' " + wellcounter + " out of " + fileList.size()/nCh) // log progress

        // get coordinates from current file
        (xyRes, width, height, coords) = getCoordinates(fileList[i])
        imageTitle = null

        chProjections = []

        // iterate over channels of a well  ********************************************************
        for (int j = 0; j < nCh; j++) {
            // set the image title (add channel sequence)
            if (imageTitle == null) {
                imageTitle = fileList[(i * nCh) + j].getName().replace(suffix, "")
                chSeq.each {
                    imageTitle += "_" + it
                }
            }

            // Project and optional BG subtraction ********************************************************
            print("\tProjecting channel " + chSeq[j] + " ... ") // log progress
            // if BF channel images
            if (j == BFch) {
                projections = openAndProject(fileList[(i * nCh) + j], "min", false)
            }
            else {
                projections = openAndProject(fileList[(i * nCh) + j], "max", bgCorrection)
            }
            chProjections.add(projections)
            println("done") // log progress
        }

        // Merge channels, if applicable
        mergedProjections = null
        channels = null
        if (chProjections.size() > 1) {
            print("\tMerging channels ") // log progress
            // merge the channels and save them in this array   ********************************************************
            mergedProjections = []
            for (int j = 0; j < chProjections[0].size(); j++) {
                chMerger = new RGBStackMerge()
                channels = new ImagePlus[nCh]
                curCal = chProjections[0][j].getCalibration().copy()
                for (int k = 0; k < nCh; k++) {
                    channels[k] = chProjections[k][j]
                }
                merge = chMerger.mergeChannels(channels, false)
                // calibrate the image
                merge.setCalibration(curCal)
                mergedProjections << merge
                print(".")
            }
            println(" done")
        }
        else {
            mergedProjections = projections
        }

        //println(IJ.currentMemory())
        // close images in channels array
        if (channels != null) {
            channels.each {
                it.close()
            }
        }



        // stitch the image     ********************************************************
        print("\tStitching ...")
        fused = stitchInMemory(fileList[i], mergedProjections, coords, imageTitle)
        println(" done")

        // change LUT of BF
        if (BFch >= 0) {
            print("\tAdjusting LUTs ... ")
            // fused is a compositeImage object
            fused.setChannelLut(LUT.createLutFromColor(Color.WHITE), BFch +1)
            //fused.updateAndDraw()
            println("done")
        }
        //fused.show()

        curCal = projections[0].getCalibration().copy()
        fused.setCalibration(curCal)


        // close projection images to free space
        projections.each {
            it.close()
        }
        IJ.run("Collect Garbage", "")
        if (!chProjections.isEmpty()) {
            chProjections.each { ch ->
                ch.each {
                    it.close()
                }
            }
        }
        IJ.run("Collect Garbage", "")
        if (mergedProjections != null) {
            mergedProjections.each {
                it.close()
            }
        }
        IJ.run("Collect Garbage", "")


        // save stitched image      ********************************************************
        String path = output.getAbsolutePath() + File.separator + fused.title + ".tif"
        println("\tSaving stitched image to: " + path)
        IJ.saveAs(fused, "Tiff", path)
        fused.close()
        IJ.run("Collect Garbage", "")


        // TODO 4 write TileConfig file?

        println("\tProcessing 'well' " + wellcounter + " done")
        wellcounter++
    } // end 'well' loop

} // end main



//  ***********         Functions      ***********//

/**
 * Create a list of projected ImagePlus objects
 * @param file = nd2 file to open
 * @param method = Projection method ("max" or "min")
 * @param subBG = Boolean do a background subtraction?
 * @return list of ImagePlus objects
 */
def openAndProject(File file, String method, Boolean subBG) {
    // open all series of container with BF **********************************
    options = new ImporterOptions()
    // import options FIXME only for one file
    options.setId(file.getAbsolutePath())
    options.setAutoscale(true)
    options.setCrop(false)
    options.setColorMode(ImporterOptions.COLOR_MODE_DEFAULT)
    options.setOpenAllSeries(true)

    // add all images to arry
    imps = BF.openImagePlus(options)

    // create maximum projection of all tiles and save them in a list
    projections = []
    imps.each {
        if (subBG) {
            curCal = it.getCalibration().copy()
            bgS = new BackgroundSubtracter()
            stack = it.getImageStack()
            for (int l = 1; l <= stack.getSize(); l++) {
                bgS.rollingBallBackground(stack.getProcessor(l), 50.0, false, false, false, false, true)
            }
            impBGS = new ImagePlus(it.title, stack)
            impBGS.setCalibration(curCal)
            projections << (ZProjector.run(impBGS, method))
        }
        else {
            projections << (ZProjector.run(it, method))
        }
        it.close()
    }
    return projections
}

/**
 * Takes care of stitching a list of imagesPlus objects
 * @param projections = imagePlus list
 * @param coords = Map as outputted by getCoordinates (-> patches)
 * @param imageTitle = String title for the image
 * @return = fused imagePlus object (not calibrated)
 */
def stitchInMemory(File file, List projections, Object coords, String imageTitle) {
    // prepare the ImageCollectionElement list for stitching
    tiles = []
    projections.eachWithIndex { mip, mipIndex ->
        ice = new ImageCollectionElement(file, mipIndex)
        ice.setImagePlus(mip)
        ice.setModel(new TranslationModel2D()) // always 2d for MIPs
        ice.setOffset([(coords[0].elements[mipIndex].x), (coords[0].elements[mipIndex].y)] as float[])
        tiles << ice
    }

    // print tile positions
    /*
    println("ice tiles:")
    tiles.each {
        println(it.getIndex() + " " + it.getOffset())
    }
    */

    // stitching (from Jan's script)
    params = new StitchingParameters()
    params.dimensionality = 2
    params.cpuMemChoice = 1 // More RAM, less time FIXME offer choice?
    params.channel1 = 0
    params.channel2 = 0
    params.timeSelect = 0
    params.checkPeaks = 5
    params.regThreshold = 0.3
    params.computeOverlap = true
    params.displayFusion = false
    params.subpixelAccuracy = true
    params.fusionMethod = 0 // Linear Blending
    params.outputDirectory = null

    stitchedList = CollectionStitchingImgLib.stitchCollection(tiles, params)


    // {stuff missing here from jans script}

    // get lists of images and models
    // FIXME: I believe it could all be done avoiding the ICE?! but what about the offset? is it calculate in the model??
    // FIXME continued: probably stitchedList converted something in the model?
    images = []
    models = []
    stitchedList.each {
        images << it.getImagePlus()
        models << it.getModel()
    }

    // start the fusion (stitching)
    noOverlap = false
    type = ImageJFunctions.wrap(images[0]).firstElement()
    fused = Fusion.fuse(type, images, models, params.dimensionality, params.subpixelAccuracy, params.fusionMethod, params.outputDirectory, noOverlap, false, params.displayFusion)

    fused.setTitle(imageTitle)

    // close images
    images.each {
        it.close()
    }
    return fused
}


/**
 * gets the coordinates from multipoint nd2 file
 * and adjust them to relative point of origin
 * inverts also te stage coordinate, only X positions (specific for Ti2 nd2?!)
 *
 * @param file = File (nd2)
 * @return pixelSize = µm/pixel
 * @return width, height = FOV dimensions in µm
 * @return patches = x/y positions (dictionary) in "reference frame" units (should be in µm)
 */
def getCoordinates(File file) {
    reader = new ImageReader()
    reader.setId(file.getAbsolutePath())
    sCount = reader.seriesCount
    pixelSize = reader.getMetadataValue("dCalibration")

    width = reader.getSizeX() * pixelSize // width in µm
    height = reader.getSizeY() * pixelSize // height in µm

    coords = [:]
    for (int i = 0; i < sCount; i++) {
        reader.setSeries(i)
        metadata = reader.getSeriesMetadata()
        xy = [x: metadata.get("X position").getAt("value"), y: metadata.get("Y position").getAt("value")]
        xy.x /= pixelSize
        xy.y /= pixelSize
        coords[i] = xy // x/y position as Double converted to pixel values
    }

    // adjust coordinates to relative 0/0 position and invert x coordinates (not y)
    xAbs = coords[0].x
    yAbs = coords[0].y
    coords.each { key, pos ->
        pos.x = -(pos.x - xAbs) // -width
        pos.y = (pos.y - yAbs) // -height
    }

    // convert width and height to pixel units for next step
    wPx = width / pixelSize
    hPx = height / pixelSize

    /*
        from Jan Eglinger's script (thanks to him)
        Takes the coordinats and creates new array with element = coords,
        and the corner points of the stage (min/max of all positions)
     */
    patches = []
    coords.each { key, pos ->
        for (patch in patches) {
            if ((pos.x >= patch.xmin - wPx) &&
                    (pos.x <= patch.xmax + wPx) &&
                    (pos.y >= patch.ymin - hPx) &&
                    (pos.y <= patch.ymax + hPx)) {
                patch.elements[key] = pos
                patch.xmin = Math.min(patch.xmin, pos.x)
                patch.xmax = Math.max(patch.xmax, pos.x)
                patch.ymin = Math.min(patch.ymin, pos.y)
                patch.ymax = Math.max(patch.ymax, pos.y)
                return // continue to next position
            }
        }
        patches << [elements: [(key):pos], xmin: pos.x, xmax: pos.x, ymin: pos.y, ymax: pos.y]
    }

    /*
    // check array
    patches[0].elements.each { key, value ->
        println("index = " + key)
        println(value.x)
    }
    */

    reader.close()
    return [pixelSize, width, height, patches]
}


/**
 * Get information about channel names from user
 * will automatically fill in information for standard channel names:
 * BF, DAPI, 470, 555, 640
 * @param fileList = list of files
 * @param nCh = number of channels
 * @return = string list of channels
 */
def getChannelIdentifier(fileList, nCh) {
    // ask user about the channel identification
    gd = new GenericDialogPlus("Channel identification")
    gd.addMessage("Channel Information: e.g. Ch1 = 555 or DAPI")
    for (int i = 0; i < nCh; i++) {
        gd.addMessage("Ch" + (i+1).toString() +": " +fileList[i].name)
    }

    // try to identify channels
    possCh = new String[nCh]
    for (int i = 0; i < nCh; i++) {
        possCh[i] = "dapi"
    }
    for (int i = 0; i < nCh; i++) {
        if (fileList[i].name.contains("DAPI")) possCh[i] = "DAPI"
        if (fileList[i].name.contains("BF")) possCh[i] = "BF"
        if (fileList[i].name.contains("470")) possCh[i] = "470"
        if (fileList[i].name.contains("555")) possCh[i] = "555"
        if (fileList[i].name.contains("640")) possCh[i] = "640"
    }

    gd.addMessage("\n")
    for (int i = 0; i < nCh; i++) {
        gd.addStringField("Channel "+(i+1).toString()+" identifier:", possCh[i])
    }
    // identify the BF channel
    if (possCh.findIndexOf {it in "BF"} != -1) {
        gd.addNumericField("BrightField Channel Number (leave if no BF)", possCh.findIndexOf {it in "BF"} + 1, 0)
    }
    else {
        gd.addNumericField("BrightField Channel Number (leave if no BF)", 0, 0)
    }
    gd.hideCancelButton()
    gd.showDialog()

    chStrings = []
    for (int i = 0; i < nCh; i++) {
        chStrings.add(gd.getNextString())
    }
    //if (chStrings.contains(gd.getNextString())) println("BF is there")
    BFch = (gd.getNextNumber() - 1) as Integer

    return [chStrings, BFch]
}