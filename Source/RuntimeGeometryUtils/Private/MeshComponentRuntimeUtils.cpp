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


RUNTIMEGEOMETRYUTILS_API void RTGUtils::UpdateDynamicMeshFromStaticMesh(UStaticMesh* StaticMesh, FDynamicMesh3& OutMesh, bool bUseMeshDescriptorsInEditor /*= false*/)
{

#if WITH_EDITORONLY_DATA
	if (bUseMeshDescriptorsInEditor)
	{
		//Potentially faster path in editor
		FMeshDescriptionToDynamicMesh Converter;
		Converter.bPrintDebugMessages = true;
		//Converter.bEnableOutputGroups = false;

		FMeshDescription* Description = StaticMesh->GetMeshDescription(0);
		//Description->DeletePolygonGroup(FPolygonGroupID(0));

		Converter.Convert(Description, OutMesh);
		return;
	}
	else
	{
		if (StaticMesh && !StaticMesh->bAllowCPUAccess)
		{
			StaticMesh->bAllowCPUAccess = true;
			UE_LOG(LogTemp, Warning,
				TEXT("RTGUtils::UpdateDynamicMeshFromStaticMesh: Overriding bAllowCPUAccess for %s in editor mode only! Update setting for asset or it will not work in a packaged game."),
				*StaticMesh->GetFName().ToString());
		}
	}
#endif

	if (!(StaticMesh && StaticMesh->bAllowCPUAccess && StaticMesh->GetRenderData() != nullptr))
	{
		UE_LOG(LogTemp, Warning, TEXT("RTGUtils::UpdateDynamicMeshFromStaticMesh: Cannot convert static mesh to FDynamicMesh3. Static mesh requires AllowCPUAccess set to true."));
		return;
	}

	OutMesh.EnableAttributes();
	FDynamicMeshNormalOverlay* Normals = OutMesh.Attributes()->PrimaryNormals();
	FDynamicMeshUVOverlay* UVs = OutMesh.Attributes()->PrimaryUV();


	if (StaticMesh->GetRenderData()->LODResources.IsValidIndex(0))
	{
		int32 NumSections = StaticMesh->GetNumSections(0);
		const FStaticMeshLODResources& LOD = StaticMesh->GetRenderData()->LODResources[0];
		const FStaticMeshVertexBuffers& VertexBuffers = LOD.VertexBuffers;

		//Copy sections like a procedural mesh method
		TArray<int32> Triangles;
		for (int32 SectionIndex = 0; SectionIndex < NumSections; SectionIndex++)
		{
			if (LOD.Sections.IsValidIndex(SectionIndex))
			{
				TArray<int32> MeshToSectionVertMap;

				const FStaticMeshSection& Section = LOD.Sections[SectionIndex];
				const uint32 OnePastLastIndex = Section.FirstIndex + Section.NumTriangles * 3;
				FIndexArrayView Indices = LOD.IndexBuffer.GetArrayView();
				for (uint32 i = Section.FirstIndex; i < OnePastLastIndex; i++)
				{
					uint32 MeshVertIndex = Indices[i];

					int32 NewIndex = MeshToSectionVertMap.Find(MeshVertIndex);
					if (NewIndex == INDEX_NONE)
					{
						OutMesh.AppendVertex(VertexBuffers.PositionVertexBuffer.VertexPosition(MeshVertIndex));
						MeshToSectionVertMap.Add(MeshVertIndex);
					}

					Normals->AppendElement(
						FVector(VertexBuffers.StaticMeshVertexBuffer.VertexTangentZ(MeshVertIndex)));
					UVs->AppendElement(VertexBuffers.StaticMeshVertexBuffer.GetVertexUV(MeshVertIndex, 0));

					Triangles.Add(MeshVertIndex);
				}
			}
		}

		//Same as Obj reader method
		for (int Index = 0; Index < Triangles.Num(); Index += 3)
		{
			auto TriangleId = OutMesh.AppendTriangle(Triangles[Index + 0], Triangles[Index + 1], Triangles[Index + 2]);

			if (Normals && Normals->IsElement(Index + 0) && Normals->IsElement(Index + 1) && Normals->IsElement(
				Index + 2))
			{
				Normals->SetTriangle(TriangleId, FIndex3i(Index + 0, Index + 1, Index + 2));
			}
			if (UVs && UVs->IsElement(Index + 0) && UVs->IsElement(Index + 1) && UVs->IsElement(Index + 2))
			{
				UVs->SetTriangle(TriangleId, FIndex3i(Index + 0, Index + 1, Index + 2));
			}
		}
	}
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
