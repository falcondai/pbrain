import math
import vtk
from events import EventHandler

class VTKSurface(vtk.vtkActor):
    """
    CLASS: VTKSurface
    DESCR: Handles a .vtk structured points (volume) or polydata (surface) file.
    """

    def set_matrix(self, registration_mat):
        print "VTKSurface.set_matrix(", registration_mat, ")!!"
        def vtkmatrix4x4_to_array(vtkmat):
            scipy_array = zeros((4,4), 'd')
            for i in range(0,4):
                for j in range(0,4):
                    scipy_array[i][j] = mat.GetElement(i,j)
            return scipy_array 

        def array_to_vtkmatrix4x4(scipy_array):
            mat = vtk.vtkMatrix4x4()
            for i in range(0,4):
                for j in range(0,4):
                    mat.SetElement(i,j, scipy_array[i][j])
            return mat

        #print "calling SetUserMatrix(", array_to_vtkmatrix4x4(registration_mat) , ")"
        mat = array_to_vtkmatrix4x4(registration_mat)
        mat.Modified()

        mat2xform = vtk.vtkMatrixToLinearTransform()
        mat2xform.SetInput(mat)
        
        print "calling SetUserTransform(", mat2xform, ")"
        self.SetUserTransform(mat2xform) # see vtk Prop3d docs
        self.Modified()
        # how do we like update the render tree or somethin..
        self.renderer.Render()

    def __init__(self, filename, renderer):

        self.renderer = renderer
        
        reader = vtk.vtkDataSetReader()
        #reader.SetFileName('/home/mcc/src/devel/extract_mri_slices/braintest2.vtk')
        reader.SetFileName(filename)

        strip = None
        if reader.IsFileStructuredPoints():
            print 'Structured points input. Generating surface from volume data.'
            strip = self.mesh_from_vol(reader.GetStructuredPointsOutput())
        elif reader.IsFilePolyData():
            print 'Poly data input.'
            strip = self.mesh_from_surf(reader.GetPolyDataOutput())
        else:
            print 'Unsupported file format. Only polydata or structured points are supported.'
        
        if not strip:
            return


        # we want to move this from its (.87 .92 .43) esque position to something more like 'the center'
        # how to do this?!?

        # ALTERNATIVELY: we want to use vtkInteractorStyleTrackballActor
        # somewhere instead of the interactor controlling the main window and 3 planes
        
        '''
        imagedata = reader.GetOutput()

        #reader.SetFileName(filename)
        cf = vtk.vtkContourFilter()
        cf.SetInput(imagedata)
        # ??? 
        cf.SetValue(0, 1)

        deci = vtk.vtkDecimatePro()
        deci.SetInput(cf.GetOutput())
        deci.SetTargetReduction(.1)
        deci.PreserveTopologyOn()


        smoother = vtk.vtkSmoothPolyDataFilter()
        smoother.SetInput(deci.GetOutput())
        smoother.SetNumberOfIterations(100)



        # XXX try to call SetScale directly on actor..
        #self.scaleTransform = vtk.vtkTransform()
        #self.scaleTransform.Identity()
        #self.scaleTransform.Scale(.1, .1, .1)
        


        #transformFilter = vtk.vtkTransformPolyDataFilter()
        #transformFilter.SetTransform(self.scaleTransform)
        #transformFilter.SetInput(smoother.GetOutput())


        #cf.SetValue(1, 2)
        #cf.SetValue(2, 3)
        #cf.GenerateValues(0, -1.0, 1.0)
        
        #deci = vtk.vtkDecimatePro()
        #deci.SetInput(cf.GetOutput())
        #deci.SetTargetReduction(0.8) # decimate_value

        normals = vtk.vtkPolyDataNormals()
        #normals.SetInput(transformFilter.GetOutput())
        normals.SetInput(smoother.GetOutput())
        normals.FlipNormalsOn()
        '''

        """
        tags = vtk.vtkFloatArray()
        tags.InsertNextValue(1.0)
        tags.InsertNextValue(0.5)
        tags.InsertNextValue(0.7)
        tags.SetName("tag")
        """

        lut = vtk.vtkLookupTable()
        lut.SetHueRange(0, 0)
        lut.SetSaturationRange(0, 0)
        lut.SetValueRange(0.2, 0.55)
        
        contourMapper = vtk.vtkPolyDataMapper()
        contourMapper.SetInputConnection(strip.GetOutputPort())
        contourMapper.SetLookupTable(lut)

        ###contourMapper.SetColorModeToMapScalars()
        ###contourMapper.SelectColorArray("tag")
        
        self.contours = vtk.vtkActor()
        self.contours.SetMapper(contourMapper)
        #if (do_wireframe):
        #self.contours.GetProperty().SetRepresentationToWireframe()
        #elif (do_surface):
        self.contours.GetProperty().SetRepresentationToSurface()
        self.contours.GetProperty().SetInterpolationToPhong()
        self.contours.GetProperty().SetOpacity(1.0)
        self.contours.GetProperty().SetAmbient(0.2)
        self.contours.GetProperty().SetDiffuse(0.5)
        self.contours.GetProperty().SetSpecular(0.1)
        self.contours.GetProperty().SetSpecularPower(0.1)

        # XXX arbitrarily setting scale to this
        #self.contours.SetScale(.1, .1,.1)

        renderer.AddActor(self.contours)
        # XXX: mcc will this work?!?

        print "PlaneWidgetsXYZ.set_image_data: setting EventHandler.set_vtkactor(self.contours)!"
        EventHandler().set_vtkactor(self.contours)

        #writer = vtk.vtkSTLWriter()
        #writer.SetFileTypeToBinary()
        #writer.SetFileName('/home/mcc/src/devel/extract_mri_slices/braintest2.stl')
        #writer.SetInput(normals.GetOutput())
        #writer.Write()
        ######################################################################
        ######################################################################
        ######################################################################

    def mesh_from_vol(self, structured_points):
        cf = vtk.vtkContourFilter()
        cf.SetInput(structured_points)
        cf.SetValue(0, 1)
        deci = vtk.vtkDecimatePro()
        deci.SetInput(cf.GetOutput())
        deci.SetTargetReduction(.1)
        deci.PreserveTopologyOn()

        smoother = vtk.vtkSmoothPolyDataFilter()
        smoother.SetInput(deci.GetOutput())
        smoother.SetNumberOfIterations(100)

        normals = vtk.vtkPolyDataNormals()
        normals.SetInput(smoother.GetOutput())
        #normals.FlipNormalsOn()
        normals.SetFeatureAngle(60.0)

        stripper = vtk.vtkStripper()
        stripper.SetInputConnection(normals.GetOutputPort())

        return stripper

    def mesh_from_surf(self, poly_data):
        normals = vtk.vtkPolyDataNormals()
        normals.SetInput(poly_data)
        normals.SetFeatureAngle(60.0)
        stripper = vtk.vtkStripper()
        stripper.SetInputConnection(normals.GetOutputPort())

        return stripper
