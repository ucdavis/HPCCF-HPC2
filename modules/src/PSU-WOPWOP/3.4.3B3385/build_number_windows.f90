function GetBuildNumber () result (build)
  character(len=128) :: build
  build="[$Revision: 3356 $]"
end function GetBuildNumber