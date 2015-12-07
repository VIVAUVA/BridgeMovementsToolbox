"""
Name:      BridgeMovements.pyt
Purpose:
Author:
Created:
Copyright: (c) 2015

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

# Import arcpy
import arcpy as ap
import numpy as np
import os
import calendar
from datetime import datetime
from datetime import timedelta

class Toolbox(object):
    def __init__(self):
        self.label = "bridge movement dev"
        self.alias = "bridgeMov_dev"
        self.tools = [BridgeMovements]


class BridgeMovements(object):

    def __init__(self):
        self.label = "Bridge movements"
        self.description = "Bridge movements"
        self.canRunInBackground = False
        return

    def getParameterInfo(self):
        # Parameter: the feature layer containing the SqueeSAR data
        sqsar_layer = ap.Parameter(
            displayName="Layer containing SqueeSAR data",
            name="sqsar_layer",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")
        sqsar_layer.value = "SqueeSAR_Demo_Data"

        # Parameter: the feature layer containing the bridges as features
        bridge_layer = ap.Parameter(
            displayName="Layer containing bridges as features",
            name="bridge_layer",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")
        bridge_layer.value = "Bridges_Demo_Data"

        # Parameter: the symbology layer containing color symbology for output
        # symbol_layer = ap.Parameter(
            # displayName="Layer containing color symbology",
            # name="symbol_layer",
            # datatype="GPFeatureLayer",
            # parameterType="Required",
            # direction="Input")
        # symbol_layer.value = "BridgeMovementsSymbology"

        # Parameter: threshold for warning situation
        warn_thresh = ap.Parameter(
            displayName="Displacement velocity to use as warning threshold [inch/year]",
            name="warn_thresh",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        warn_thresh.value = 1

        # Parameter: threshold for critical situation
        critical_thresh = ap.Parameter(
            displayName="Displacement velocity to use as critical threshold [inch/month]",
            name="critical_thresh",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        critical_thresh.value = .5

        # Parameter: output layer that will contain the results of the bridge
        # analysis
        out_layer = ap.Parameter(
            displayName="Output layer showing the results",
            name="out_layer",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Output")

        # params = [sqsar_layer, bridge_layer, symbol_layer, warn_thresh, critical_thresh, out_layer]
        params = [sqsar_layer, bridge_layer, warn_thresh, critical_thresh, out_layer]

        return params

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        return

    def udateMessages(self, parameters, messages):
        return

    def execute(self, parameters, messages):
        # The approach:
        # - define a search cursor to loop through the bridges
        # - extract the SqueeSAR points falling within the area identify by
        #   each bridge
        # - analyze these points

        # Tell ArcGIS that he/she can overwrite outputs
        ap.env.overwriteOutput = True

        # INPUTS
        # Get the SqueeSAR layer
        sqsar_layer = parameters[0].valueAsText
        # ap.AddMessage("SqueeSAR layer: {0}".format(sqsar_layer))

        # Get the bridge layer
        bridge_layer = parameters[1].valueAsText
        # ap.AddMessage("Bridge layer: {0}".format(bridge_layer))
        bridge_x = ap.da.FeatureClassToNumPyArray(bridge_layer, ["SHAPE@X"])
        bridge_y = ap.da.FeatureClassToNumPyArray(bridge_layer, ["SHAPE@Y"])
        bridge_x = bridge_x.astype(np.float)
        bridge_y = bridge_y.astype(np.float)

        # Get the symbology layer
        # symbol_layer = parameters[2].valueAsText
        # ap.AddMessage("Symbology layer: {0}".format(symbol_layer))

        # Get the warning threshold
        warn_thresh = parameters[2].valueAsText
        warn_thresh = float(warn_thresh)
        ap.AddMessage("Warning threshold: {0} inches per year".format(warn_thresh))

        # Get the critical threshold
        critical_thresh = parameters[3].valueAsText
        critical_thresh = float(critical_thresh)
        ap.AddMessage("Critical threshold: {0} inches per month".format(critical_thresh))

        # OUTPUT
        out_layer = parameters[4].valueAsText
        ap.AddMessage("Output layer: {0}".format(out_layer))

        # Copy the bridge layer to the output layer
        try:
            ap.management.MakeFeatureLayer(bridge_layer, out_layer)
        except:
            ap.AddError("Failure to copy {0} to {1}".format(bridge_layer, out_layer))
            return

        STAT_FLD = u"SAR_STAT" # Define the status field
        COND_FLD = u"SAR_COND" # Define text status field
        X_FLD = u"COOR_X"      # Define x-coordinate field for bridge
        Y_FLD = u"COOR_Y"      # Define y-coordinate field for bridge

        BRDG = u"BRDG"   # Define on-bridge field for squeesar
        defaultValue = str(0)      # Not on bridge

        # Check if the status field already exists in the output layer
        if ap.ListFields(out_layer, STAT_FLD):
            ap.AddWarning("Overwriting existing {0} field values".format(STAT_FLD))
        else:
            # Add the status field to the shape file
            try:
                ap.management.AddField(out_layer, STAT_FLD, "FLOAT")
            except:
                ap.AddError("Failure to create {0} field in {1}".format(STAT_FLD, out_layer))
                return

        # Check if the condition field already exists in the output layer
        if ap.ListFields(out_layer, COND_FLD):
            ap.AddWarning("Overwriting existing {0} field values".format(COND_FLD))
        else:
            # Add the condition field to the shape file
            try:
                ap.management.AddField(out_layer, COND_FLD, "TEXT")
            except:
                ap.AddError("Failure to create {0} field in {1}".format(COND_FLD, out_layer))
                return

        # Check if the x-coordinate field already exists in the output layer
        if ap.ListFields(out_layer, X_FLD):
            ap.AddWarning("Overwriting existing {0} field values".format(X_FLD))
        else:
            # Add the x-coordinate field to the shape file
            try:
                ap.management.AddField(out_layer, X_FLD, "FLOAT")
            except:
                ap.AddError("Failure to create {0} field in {1}".format(X_FLD, out_layer))
                return

        # Check if the y-coordinate field already exists in the output layer
        if ap.ListFields(out_layer, Y_FLD):
            ap.AddWarning("Overwriting existing {0} field values".format(Y_FLD))
        else:
            # Add the y-coordinate field to the shape file
            try:
                ap.management.AddField(out_layer, Y_FLD, "FLOAT")
            except:
                ap.AddError("Failure to create {0} field in {1}".format(Y_FLD, out_layer))
                return

        # Check if the on-bridge field already exists in the squeesar layer
        if ap.ListFields(sqsar_layer, BRDG):
            ap.AddWarning("Overwriting existing {0} field values".format(BRDG))
        else:
            # Add the on-bridge field to the shape file
            try:
                ap.management.AddField(sqsar_layer, BRDG, "TEXT")
                ap.AssignDefaultToField_management(sqsar_layer, BRDG, defaultValue)
            except Exception as e:
                ap.AddError("Failure to create {0} field in {1}: {2}".format(BRDG, sqsar_layer, e))
                return

        # Get the list of fields with date format in the SqueeSAR layer
        # NOTE: length of date fields will give number of readings
        date_fields = [f.name for f in ap.ListFields(sqsar_layer, u"D*")]
        # ap.AddMessage("date_fields: {0}".format(date_fields))

        # Create an update cursor to navigate trhough the bridges (features)
        bridge_cursor = ap.UpdateCursor(out_layer)
        # initialize blank string for bad bridge id
        bad_id = ''

        # For each bridge in the layer
        b = bridge_cursor.next()
        bridge_index = 0 # to index through array of bridge coordinates
        while b:
            # Get bridge ID
            b_id = str(b.getValue("FID"))

            # Add message
            ap.AddMessage("--- Processing bridge {0} ---".format(b_id))

            # Get bridge coordinates
            coord_x = bridge_x[bridge_index]
            coord_y = bridge_y[bridge_index]
            bridge_index += 1
            ap.AddMessage("This bridge is located at ({0},{1}) feet".format(coord_x,coord_y))

            # Select the current bridge
            ap.management.SelectLayerByAttribute(out_layer,
                                                 "NEW_SELECTION",
                                                 "FID = " + b_id)

            # Add squeesar points within the bridge to the selection
            ap.management.SelectLayerByLocation(sqsar_layer,
                                                "WITHIN",
                                                out_layer)

            # Remove the bridge from the selection. This is needed because I
            # couldn't find a way to replace the selection with the previous
            # command, non even using "NEW_SELECTION"
            ap.management.SelectLayerByAttribute(out_layer,
                                                 "REMOVE_FROM_SELECTION",
                                                 "FID = " + b_id)
            # initialize status
            status = 0 # good bridge
            status_strng = "GOOD"

            # Get all heights to check if points are on the bridge or not
            ssar_height = ap.da.FeatureClassToNumPyArray(sqsar_layer,['HEIGHT'])
            ssar_height = ssar_height.astype(np.float)
            # ap.AddMessage("Number of points before threshold: {0}".format(ssar_height.size))

            # Check if selection is empty
            if ssar_height.size == 0: # no SQUEE-SAR data within bridge
                # Update the bridge status field value
                b.setValue(STAT_FLD, status)
                b.setValue(COND_FLD, status_strng)
                b.setValue(X_FLD, coord_x)
                b.setValue(Y_FLD, coord_y)
                bridge_cursor.updateRow(b)

                # Next bridge
                b = bridge_cursor.next()
                continue

            # Remove squeesar points that are not on the bridge
            # Check to see if any points are more than 6 meters below highest
            # point - assuming highest point is on the bridge

            height_median = np.median(ssar_height)
            # ap.AddMessage("Median height is: {0}".format(height_median))
            height_threshold = height_median - 3 # anything over 3 meters below median point may be on road
            # ap.AddMessage("Threshold height is: {0}".format(height_threshold))

### SHOULD ALSO CHECK IF IT IS ABOVE THE HEIGHT (+-3m)

            ap.management.SelectLayerByAttribute(sqsar_layer,'SUBSET_SELECTION',
                                               '"HEIGHT" >= ' + str(height_threshold))

            # http://gis.stackexchange.com/questions/88320/using-variable-to-calculate-field-with-arcpy-as-script-in-toolbox
            defaultValue = str(1)
            # these points are on the bridge, so set them to value of 1
            ap.CalculateField_management(sqsar_layer, BRDG, "'" + defaultValue + "'", "PYTHON")

            # Get all velocities
            ssar_vel = ap.da.FeatureClassToNumPyArray(sqsar_layer,['VEL'])
            ssar_vel = ssar_vel.astype(np.float)
            ap.AddMessage("Number of points on bridge: {0}".format(ssar_vel.size))

            # Check again if selection is empty
            if ssar_vel.size == 0: # no SQUEE-SAR data within bridge
                # Update the bridge status field value
                b.setValue(STAT_FLD, status)
                b.setValue(COND_FLD, status_strng)
                b.setValue(X_FLD, coord_x)
                b.setValue(Y_FLD, coord_y)
                bridge_cursor.updateRow(b)

                # Next bridge
                b = bridge_cursor.next()
                bridge_index += 1
                continue

            # Get all displacements
            ssar_disp = ap.da.FeatureClassToNumPyArray(sqsar_layer,date_fields)

            # The following is code to check how far apart the data is spaced.
            # We want to take points that are about a month apart, in order to
            # calculate the displacement/month.

            # http://code.activestate.com/recipes/577274-subtract-or-add-a-month-to-a-datetimedate-or-datet/
            # add_month method by: Garel Alex
            def add_month(date):
                """add one month to date, maybe falling to last day of month

                :param datetime.datetime date: the date

                ::
                 >>> add_month(datetime(2014,1,31))
                datetime.datetime(2014, 2, 28, 0, 0)
                 >>> add_month(datetime(2014,12,30))
                datetime.datetime(2015, 1, 30, 0, 0)
                 """
                # number of days this month
                month_days = calendar.monthrange(date.year, date.month)[1]
                candidate = date + timedelta(days=month_days)
                # but maybe we are a month too far
                if candidate.day != date.day:
                        # go to last day of next month,
                        # by getting one day before begin of candidate month
                        return candidate.replace(day=1) - timedelta(days=1)
                else:
                        return candidate

            start_date = date_fields[0]
            year = int(start_date[1:5])
            month = int(start_date[5:7])
            day = int(start_date[7:9])
            t0 = datetime(year,month,day)
            t1 = add_month(t0)

            # initialize list of average velocities - this is a list of lists
            vel_list = []

### WHAT IF THERE ARE MORE THAN 61 DATES? SHOULD THIS BE [1:]?

            for date in date_fields[1:61]:

                year = int(date[1:5])
                month = int(date[5:7])
                day = int(date[7:9])
                t = datetime(year,month,day)

                if t < t1:
                        continue # time of this entry is not a month away from start date
                else:
                        months_apart = ((t - t0).days)/30 # to convert into months
                        vel = np.absolute((ssar_disp[date] - ssar_disp[start_date])/months_apart)
                        vel_list.append(vel)

                        start_date = date # new start date
                        t0 = t
                        t1 = add_month(t0)

### WHAT ABOUT DOING A LINEAR INTERPOLATION OF THE POINTS BETWEEN T AND T0?
### IT WOULD REDUCE THE NOISE IF THERE ARE MORE THAN 2 POINTS.

                        if np.max(vel) >= 25.4*critical_thresh:
                                status = 2.5
                                status_strng = "CRITICAL"
                                ap.AddMessage("STATUS: CRITICAL")
                                bad_id = b_id
                                break
                        else:
                                continue

            if status < 2.5:
                # ap.AddMessage("List of average velocities: {0}".format(vel_list[0]))
                vel_array = np.empty([ssar_vel.size, len(vel_list)])
                for j in range(0,len(vel_list)):
                        vel_array[0:ssar_vel.size,j] = vel_list[j] # put all average velocities in one numpy array
                # ap.AddMessage("Array of average velocities: {0}".format(vel_array))

                # if change over 1 year in any disp is over warning threshold (default 1"/year), raise warning
                # must convert to mm
                vel_medians = np.median(vel_array, axis=1) # calculate medians across time series of monthly velocities
                vel_medians = vel_medians*12 	       # convert to mm per year
                if np.amax(vel_medians) >= 25.4*warn_thresh:
                        status = 1.5 # bridge indicates warning
                        status_strng = "WARNING"
                        ap.AddMessage("STATUS: Warning")

                # If status of bridge is good, add message here stating so
                if status == 0:
                        ap.AddMessage("STATUS: Good")

            # Update the bridge status field value
            b.setValue(STAT_FLD, status)
            b.setValue(COND_FLD, status_strng)
            b.setValue(X_FLD, coord_x)
            b.setValue(Y_FLD, coord_y)
            bridge_cursor.updateRow(b)

            # Next bridge
            b = bridge_cursor.next()

        # Clear cursor
        del b
        del bridge_cursor

        # Define the location the file layer containing the graduated values
        # symbology
        # GV_FIL = os.path.join(os.path.dirname(__file__),
                              # "graduated_colors_polyg.lyr")
        # ap.AddMessage("Using {0} as symbology file".format(GV_FIL))

        # Modify the symbology of the new layer to identify the top percentile
        # ap.management.SelectLayerByLocation(GV_FIL,
                                                # "WITHIN",
                                                # out_layer)

        # Symbology
        symbol_layer = os.path.dirname(os.path.realpath(__file__)) +u"\BridgeMovementsSymbology.lyr"

        sym_layer_obj = ap.mapping.Layer(symbol_layer)
        # Fixing broken link of the symbology layer
        new_path = os.path.dirname(os.path.realpath(__file__)) + "\\BridgeMovementsDemoData\\"
        sym_layer_obj.replaceDataSource(new_path,"SHAPEFILE_WORKSPACE", "Bridges_Demo_Data")
        ap.management.ApplySymbologyFromLayer(out_layer, sym_layer_obj)
        out_layer_obj = ap.mapping.Layer(out_layer)
        if out_layer_obj.symbologyType == "GRADUATED_COLORS":
                out_layer_obj.symbology.valueField = STAT_FLD
                out_layer_obj.symbology.classBreakValues = [0.0, 1.0, 2.0, 3.0]
                out_layer_obj.symbology.classBreakLabels = ["OK", "WARNING", "CRITICAL"]

        # clearing last selection
        ap.management.SelectLayerByAttribute(sqsar_layer,'REMOVE_FROM_SELECTION',
                                               '"HEIGHT" >= ' + str(height_threshold))
        # Placing labels on the data
        # From https://gist.github.com/avaccari/80c34d561d013483bec9
        # ApLabels.py by Andrea Vaccari

        # Load the reference to the current map document.
        # Caveat: not sure if it will work if the toolbox is running in background.
        # This is because background execution is performed within a separate kernel.
        mxd = ap.mapping.MapDocument('CURRENT')

        # Get list of dataframes. There might be more than one dataframe in a map document. You
        # can use the wildcard to specify which one is of interest. Here we are "blindly" extracting
        # the first element in the returned list
        dfs = ap.mapping.ListDataFrames(mxd, "")
        df = dfs[0]

        # Get list of layers. Once again there might be more than one layer within the selected
        # dataframe. Wildcard can be used to specify the subset or specific one. We are picking
        # the layer that holds our squeesar data, and the labels for that data
        lyrs = ap.mapping.ListLayers(mxd,sqsar_layer, df) # (map_document_or_layer, {wildcard}, {data_frame})
        lyr = lyrs[0]

        # Check if the layer supports label classes
        if lyr.supports("LABELCLASSES"):
                # Get list of labelClasses. There might be more than one label class.
                # We are selecting/defining the first.
                lcs = lyr.labelClasses
                lc = lcs[0]

                lc.SQLQuery = '"BRDG" > ' + str(0)

                # Define the expression to use. Here we are using the attributes [HEIGHT]
                # and [VEL] from the display data table.
                # Notes:
                # - switching ' with " didn't work
                # - inserting escaped characters such as '\n' did not work
                lc.expression = '[HEIGHT] + "m, " + [VEL] + "mm/yr"'

                # Turn the labels on for the entire layer.
                lyr.showLabels = True

                # Refresh the current view
                ap.RefreshActiveView()

        # zooming to bad bridge!
        if bad_id:
                ap.AddMessage("Zooming to last identified bad bridge")
                ap.SelectLayerByAttribute_management(bridge_layer,"NEW_SELECTION","FID = " + bad_id )
                # mxd = ap.mapping.MapDocument('CURRENT')
                df = ap.mapping.ListDataFrames(mxd, "Layers") [0]
                df.zoomToSelectedFeatures()
                ap.RefreshActiveView()

        # Clear selections
        ap.management.SelectLayerByAttribute(bridge_layer, "CLEAR_SELECTION")
        ap.management.SelectLayerByAttribute(sqsar_layer, "CLEAR_SELECTION")
        ap.management.SelectLayerByAttribute(out_layer, "CLEAR_SELECTION")

        return