# Project Stitch Merge ND2 in Fiji
Allows you to MIP and stitch multipoint nd2 files.
If you have several channels saved as individual files,
this script will merge these channels.

## File structure
A folder containing several nd2 files. 
E.g. one nd2 file containing all tiles of a multi-well plate well for a channel.

## Usage
Drag and drop the script onto Fiji and run the script.
You will be asked to:
- specify a source and destiantion folder
- specify the number of acquired channels
- perform a rolling-ball background subtraction (optional)

If more than 1 channel, a second prompt will ask you about the channels.
Make sure that the provided information matches the channel information in the file name.
Specify also the bright-field channel if acquired.

## Further information
The MIP-ing and stitching is done in memory (no temporary files).
Optional backround subtraction with 50px rolling ball is done the z-stacks.
