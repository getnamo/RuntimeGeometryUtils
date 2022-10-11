
#pragma once

#include "DynamicMesh/DynamicMeshAABBTree3.h"
#include "Components/DynamicMeshComponent.h"
#include "Interfaces/Interface_CollisionDataProvider.h"
#include "ShapeApproximation/SimpleShapeSet3.h"
#include "RuntimeDynamicMeshComponent.generated.h"

/**
 * URuntimeDynamicMeshComponent extends USimpleDynamicMeshComponent with Physics/Collision Support.
 * Complex-as-Simple and externally-generated Simple Collision are available.
 * 
 */
UCLASS(ShowCategories=(Physics, Collision))
class RUNTIMEGEOMETRYUTILS_API URuntimeDynamicMeshComponent : public UDynamicMeshComponent
{
	GENERATED_BODY()

public:
	URuntimeDynamicMeshComponent();


	UPROPERTY(EditAnywhere, Category = "Runtime Dynamic Mesh")
	bool bUseComplexAsSimpleCollision = true;

	void SetSimpleCollisionGeometry(const UE::Geometry::FSimpleShapeSet3d& SimpleShapes, bool bDeferCollisionUpdate = false);
	void SetSimpleCollisionGeometry(UE::Geometry::FSimpleShapeSet3d&& SimpleShapes, bool bDeferCollisionUpdate = false);



	//
	// USimpleDynamicMeshComponent API overrides
	//
public:
	virtual void NotifyMeshUpdated() override;

	//
	// Component Physics API overrides and IInterface_CollisionDataProvider
	//
public:
	virtual bool GetPhysicsTriMeshData(struct FTriMeshCollisionData* CollisionData, bool InUseAllTriData) override;
	virtual bool ContainsPhysicsTriMeshData(bool InUseAllTriData) const override;
	virtual bool WantsNegXTriMesh() override;
	virtual UBodySetup* GetBodySetup() override;

protected:

	//UPROPERTY(EditAnywhere, Category = "Runtime Dynamic Mesh")
	//bool bUseAsyncCooking = true;		// todo: support

	void RegenerateCollision_Immediate();
	//void RegenerateCollision_Async();   // todo: implement

	UE::Geometry::FSimpleShapeSet3d SimpleCollisionShapes;

};