##
## Auto Generated makefile by CodeLite IDE
## any manual changes will be erased      
##
## Debug
ProjectName            :=HPCcoursework
ConfigurationName      :=Debug
WorkspacePath          :=/home/abd17/Documents/hpc-exercises
ProjectPath            :=/home/abd17/Documents/hpc-exercises/HPCcoursework
IntermediateDirectory  :=./Debug
OutDir                 := $(IntermediateDirectory)
CurrentFileName        :=
CurrentFilePath        :=
CurrentFileFullPath    :=
User                   :=Ananya Dubey
Date                   :=14/03/20
CodeLitePath           :=/home/abd17/.codelite
LinkerName             :=/usr/bin/g++
SharedObjectLinkerName :=/usr/bin/g++ -shared -fPIC
ObjectSuffix           :=.o
DependSuffix           :=.o.d
PreprocessSuffix       :=.i
DebugSwitch            :=-g 
IncludeSwitch          :=-I
LibrarySwitch          :=-l
OutputSwitch           :=-o 
LibraryPathSwitch      :=-L
PreprocessorSwitch     :=-D
SourceSwitch           :=-c 
OutputFile             :=$(IntermediateDirectory)/$(ProjectName)
Preprocessors          :=
ObjectSwitch           :=-o 
ArchiveOutputSwitch    := 
PreprocessOnlySwitch   :=-E
ObjectsFileList        :="HPCcoursework.txt"
PCHCompileFlags        :=
MakeDirCommand         :=mkdir -p
LinkOptions            :=  
IncludePath            :=  $(IncludeSwitch). $(IncludeSwitch). 
IncludePCH             := 
RcIncludePath          := 
Libs                   := $(LibrarySwitch)lapack $(LibrarySwitch)blas 
ArLibs                 :=  "lapack" "blas" 
LibPath                := $(LibraryPathSwitch). 

##
## Common variables
## AR, CXX, CC, AS, CXXFLAGS and CFLAGS can be overriden using an environment variables
##
AR       := /usr/bin/ar rcu
CXX      := /usr/bin/g++
CC       := /usr/bin/gcc
CXXFLAGS :=  -g -O0 -Wall $(Preprocessors)
CFLAGS   :=  -g -O0 -Wall $(Preprocessors)
ASFLAGS  := 
AS       := /usr/bin/as


##
## User defined environment variables
##
CodeLiteDir:=/usr/share/codelite
Objects0=$(IntermediateDirectory)/test.cpp$(ObjectSuffix) $(IntermediateDirectory)/LidDrivenCavity.cpp$(ObjectSuffix) $(IntermediateDirectory)/LidDrivenCavitySolver.cpp$(ObjectSuffix) 



Objects=$(Objects0) 

##
## Main Build Targets 
##
.PHONY: all clean PreBuild PrePreBuild PostBuild MakeIntermediateDirs
all: $(OutputFile)

$(OutputFile): $(IntermediateDirectory)/.d $(Objects) 
	@$(MakeDirCommand) $(@D)
	@echo "" > $(IntermediateDirectory)/.d
	@echo $(Objects0)  > $(ObjectsFileList)
	$(LinkerName) $(OutputSwitch)$(OutputFile) @$(ObjectsFileList) $(LibPath) $(Libs) $(LinkOptions)

MakeIntermediateDirs:
	@test -d ./Debug || $(MakeDirCommand) ./Debug


$(IntermediateDirectory)/.d:
	@test -d ./Debug || $(MakeDirCommand) ./Debug

PreBuild:


##
## Objects
##
$(IntermediateDirectory)/test.cpp$(ObjectSuffix): test.cpp $(IntermediateDirectory)/test.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/abd17/Documents/hpc-exercises/HPCcoursework/test.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/test.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/test.cpp$(DependSuffix): test.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/test.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/test.cpp$(DependSuffix) -MM test.cpp

$(IntermediateDirectory)/test.cpp$(PreprocessSuffix): test.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/test.cpp$(PreprocessSuffix) test.cpp

$(IntermediateDirectory)/LidDrivenCavity.cpp$(ObjectSuffix): LidDrivenCavity.cpp $(IntermediateDirectory)/LidDrivenCavity.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/abd17/Documents/hpc-exercises/HPCcoursework/LidDrivenCavity.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/LidDrivenCavity.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/LidDrivenCavity.cpp$(DependSuffix): LidDrivenCavity.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/LidDrivenCavity.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/LidDrivenCavity.cpp$(DependSuffix) -MM LidDrivenCavity.cpp

$(IntermediateDirectory)/LidDrivenCavity.cpp$(PreprocessSuffix): LidDrivenCavity.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/LidDrivenCavity.cpp$(PreprocessSuffix) LidDrivenCavity.cpp

$(IntermediateDirectory)/LidDrivenCavitySolver.cpp$(ObjectSuffix): LidDrivenCavitySolver.cpp $(IntermediateDirectory)/LidDrivenCavitySolver.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/abd17/Documents/hpc-exercises/HPCcoursework/LidDrivenCavitySolver.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/LidDrivenCavitySolver.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/LidDrivenCavitySolver.cpp$(DependSuffix): LidDrivenCavitySolver.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/LidDrivenCavitySolver.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/LidDrivenCavitySolver.cpp$(DependSuffix) -MM LidDrivenCavitySolver.cpp

$(IntermediateDirectory)/LidDrivenCavitySolver.cpp$(PreprocessSuffix): LidDrivenCavitySolver.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/LidDrivenCavitySolver.cpp$(PreprocessSuffix) LidDrivenCavitySolver.cpp


-include $(IntermediateDirectory)/*$(DependSuffix)
##
## Clean
##
clean:
	$(RM) -r ./Debug/


