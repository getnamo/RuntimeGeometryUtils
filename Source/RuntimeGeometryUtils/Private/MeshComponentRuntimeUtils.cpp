#include "MeshComponentRuntimeUtils.h"

#include "DynamicMeshAttributeSet.h"
#include "MeshNormals.h"

#include "DynamicMeshToMeshDescription.h"
#include "MeshDescriptionToDynamicMesh.h"
#include "StaticMeshAttributes.h"
#include "KismetProceduralMeshLibrary.h"






void RTGUtils::UpdateStaticMeshFromDynamicMesh(
	UStaticMesh* StaticMesh,
	const FDynamicMesh3* Mesh)
{
	FMeshDescription MeshDescription;
	FStaticMeshAttributes StaticMeshAttributes(MeshDescription);
	StaticMeshAttributes.Register();

	FDynamicMeshToMeshDescription Converter;
	Converter.Convert(Mesh, MeshDescription);

	// todo: vertex color support

	//UStaticMesh* StaticMesh = NewObject<UStaticMesh>(Component);
	//FName MaterialSlotName = StaticMesh->AddMaterial(MyMaterial);

	// Build the static mesh render data, one FMeshDescription* per LOD.
	TArray<const FMeshDescription*> MeshDescriptionPtrs;
	MeshDescriptionPtrs.Emplace(&MeshDescription);
	StaticMesh->BuildFromMeshDescriptions(MeshDescriptionPtrs);
}


RUNTIMEGEOMETRYUTILS_API void RTGUtils::UpdateDynamicMeshFromStaticMesh(UStaticMesh* StaticMesh, FDynamicMesh3& OutMesh)
{
	FMeshDescriptionToDynamicMesh Converter;
	Converter.bPrintDebugMessages = true;
	//Converter.bEnableOutputGroups = false;

	FMeshDescription* Description = StaticMesh->GetMeshDescription(0);
	
	//Description->DeletePolygonGroup(FPolygonGroupID(0));

	Converter.Convert(Description, OutMesh);
		
	/*
	//Raw method variant
	
	int32 LODIndex = 0;
	int32 NumSections = StaticMesh->GetNumSections(LODIndex);
	for (int32 SectionIndex = 0; SectionIndex < NumSections; SectionIndex++)
	{
		// Buffers for copying geom data
		TArray<FVector> Vertices;
		TArray<int32> Triangles;
		TArray<FVector> Normals;
		TArray<FVector2D> UVs;
		Todo: support multiple UVs?
		//TArray<FVector2D> UVs1;
		//TArray<FVector2D> UVs2;
		//TArray<FVector2D> UVs3;
		TArray<FProcMeshTangent> Tangents;

		// Get geom data from static mesh
		UKismetProceduralMeshLibrary::GetSectionFromStaticMesh(StaticMesh, LODIndex, SectionIndex, Vertices, Triangles, Normals, UVs, Tangents);

		// Append vertices
		for (FVector& Vertex : Vertices)
		{
			OutMesh.AppendVertex(Vertex);
		}

		OutMesh.EnableAttributes();

		FDynamicMeshNormalOverlay* OutNormals = OutMesh.Attributes()->PrimaryNormals();
		FDynamicMeshUVOverlay* OutUVs = OutMesh.Attributes()->PrimaryUV();
		

		// Normals
		for (FVector& Normal : Normals)
		{
			OutNormals->AppendElement(Normal);
		}

		// UVs
		for (FVector2D& UV : UVs)
		{
			OutUVs->AppendElement(UV);
		}

		//Triangles
		for (int i = 0; (i+2) <Triangles.Num(); i = i+3)
		{
			FIndex3i TriangleI3;
			TriangleI3.A = Triangles[i];
			TriangleI3.B = Triangles[i + 1];
			TriangleI3.C = Triangles[i + 2];

			int32 TriangleId = OutMesh.AppendTriangle(TriangleI3);
			OutUVs->SetTriangle(TriangleId, TriangleI3);
		}
	}*/
}

void RTGUtils::UpdatePMCFromDynamicMesh_SplitTriangles(
	UProceduralMeshComponent* Component, 
	const FDynamicMesh3* Mesh,
	bool bUseFaceNormals,
	bool bInitializeUV0,
	bool bInitializePerVertexColors,
	bool bCreateCollision)
{
	Component->ClearAllMeshSections();

	int32 NumTriangles = Mesh->TriangleCount();
	int32 NumVertices = NumTriangles * 3;

	TArray<FVector> Vertices, Normals;
	Vertices.SetNumUninitialized(NumVertices);
	Normals.SetNumUninitialized(NumVertices);

	FMeshNormals PerVertexNormals(Mesh);
	bool bUsePerVertexNormals = false;
	const FDynamicMeshNormalOverlay* NormalOverlay = nullptr;
	if (Mesh->HasAttributes() == false && bUseFaceNormals == false)
	{
		PerVertexNormals.ComputeVertexNormals();
		bUsePerVertexNormals = true;
	}
	else if (Mesh->HasAttributes())
	{
		NormalOverlay = Mesh->Attributes()->PrimaryNormals();
	}

	const FDynamicMeshUVOverlay* UVOverlay = (Mesh->HasAttributes()) ? Mesh->Attributes()->PrimaryUV() : nullptr;
	TArray<FVector2D> UV0;
	if (UVOverlay && bInitializeUV0)
	{
		UV0.SetNum(NumVertices);
	}

	TArray<FLinearColor> VtxColors;
	bool bUsePerVertexColors = false;
	if (bInitializePerVertexColors && Mesh->HasVertexColors())
	{
		VtxColors.SetNum(NumVertices);
		bUsePerVertexColors = true;
	}

	TArray<FProcMeshTangent> Tangents;		// not supporting this for now

	TArray<int32> Triangles;
	Triangles.SetNumUninitialized(NumTriangles*3);

	FVector3d Position[3];
	FVector3f Normal[3];
	FVector2f UV[3];
	int32 BufferIndex = 0;
	for (int32 tid : Mesh->TriangleIndicesItr())
	{
		int32 k = 3 * (BufferIndex++);

		FIndex3i TriVerts = Mesh->GetTriangle(tid);

		Mesh->GetTriVertices(tid, Position[0], Position[1], Position[2]);
		Vertices[k] = (FVector)Position[0];
		Vertices[k+1] = (FVector)Position[1];
		Vertices[k+2] = (FVector)Position[2];


		if (bUsePerVertexNormals)
		{
			Normals[k] = (FVector)PerVertexNormals[TriVerts.A];
			Normals[k+1] = (FVector)PerVertexNormals[TriVerts.B];
			Normals[k+2] = (FVector)PerVertexNormals[TriVerts.C];
		}
		else if (NormalOverlay != nullptr && bUseFaceNormals == false)
		{
			NormalOverlay->GetTriElements(tid, Normal[0], Normal[1], Normal[2]);
			Normals[k] = (FVector)Normal[0];
			Normals[k+1] = (FVector)Normal[1];
			Normals[k+2] = (FVector)Normal[2];
		}
		else
		{
			FVector3d TriNormal = Mesh->GetTriNormal(tid);
			Normals[k] = (FVector)TriNormal;
			Normals[k+1] = (FVector)TriNormal;
			Normals[k+2] = (FVector)TriNormal;
		}

		if (UVOverlay != nullptr && UVOverlay->IsSetTriangle(tid))
		{
			UVOverlay->GetTriElements(tid, UV[0], UV[1], UV[2]);
			UV0[k] = (FVector2D)UV[0];
			UV0[k+1] = (FVector2D)UV[1];
			UV0[k+2] = (FVector2D)UV[2];
		}

		if (bUsePerVertexColors)
		{
			VtxColors[k] = (FLinearColor)Mesh->GetVertexColor(TriVerts.A);
			VtxColors[k+1] = (FLinearColor)Mesh->GetVertexColor(TriVerts.B);
			VtxColors[k+2] = (FLinearColor)Mesh->GetVertexColor(TriVerts.C);
		}

		Triangles[k] = k;
		Triangles[k+1] = k+1;
		Triangles[k+2] = k+2;
	}

	Component->CreateMeshSection_LinearColor(0, Vertices, Triangles, Normals, UV0, VtxColors, Tangents, bCreateCollision);
}
