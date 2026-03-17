#Tue Mar 17 23:24:47 2026"""Automatically generated. DO NOT EDIT please"""
from GaudiKernel.Proxy.Configurable import *

class BasicDigitizer( ConfigurableAlgorithm ) :
  __slots__ = { 
    'ExtraInputs' : set(),
    'ExtraOutputs' : set(),
    'OutputLevel' : 0,
    'Enable' : True,
    'ErrorMax' : 1,
    'AuditInitialize' : False,
    'AuditReinitialize' : False,
    'AuditRestart' : False,
    'AuditExecute' : False,
    'AuditFinalize' : False,
    'AuditStart' : False,
    'AuditStop' : False,
    'Timeline' : True,
    'MonitorService' : 'MonitorSvc',
    'RegisterForContextService' : False,
    'Cardinality' : 0,
    'NeededResources' : [  ],
    'Asynchronous' : False,
    'FilterCircularDependencies' : True,
    'InputCollection' : 'SiTargetHitsMIP',
    'OutputCollection' : 'SiTargetDigiHits',
    'Threshold' : 0.50000000,
    'DebugFrequency' : 500,
  }
  _propertyDocDct = { 
    'ExtraInputs' : """  [DataHandleHolderBase<PropertyHolder<CommonMessaging<implements<IAlgorithm,IDataHandleHolder,IProperty,IStateful> > > >] """,
    'ExtraOutputs' : """  [DataHandleHolderBase<PropertyHolder<CommonMessaging<implements<IAlgorithm,IDataHandleHolder,IProperty,IStateful> > > >] """,
    'OutputLevel' : """ output level [Gaudi::Algorithm] """,
    'Enable' : """ should the algorithm be executed or not [Gaudi::Algorithm] """,
    'ErrorMax' : """ [[deprecated]] max number of errors [Gaudi::Algorithm] """,
    'AuditInitialize' : """ trigger auditor on initialize() [Gaudi::Algorithm] """,
    'AuditReinitialize' : """ trigger auditor on reinitialize() [Gaudi::Algorithm] """,
    'AuditRestart' : """ trigger auditor on restart() [Gaudi::Algorithm] """,
    'AuditExecute' : """ trigger auditor on execute() [Gaudi::Algorithm] """,
    'AuditFinalize' : """ trigger auditor on finalize() [Gaudi::Algorithm] """,
    'AuditStart' : """ trigger auditor on start() [Gaudi::Algorithm] """,
    'AuditStop' : """ trigger auditor on stop() [Gaudi::Algorithm] """,
    'Timeline' : """ send events to TimelineSvc [Gaudi::Algorithm] """,
    'MonitorService' : """ name to use for Monitor Service [Gaudi::Algorithm] """,
    'RegisterForContextService' : """ flag to enforce the registration for Algorithm Context Service [Gaudi::Algorithm] """,
    'Cardinality' : """ how many clones to create - 0 means algo is reentrant [Gaudi::Algorithm] """,
    'NeededResources' : """ named resources needed during event looping [Gaudi::Algorithm] """,
    'Asynchronous' : """ whether algorithm is asynchronous and uses Boost Fiber to suspend while offloaded code is running. [Gaudi::Algorithm] """,
    'FilterCircularDependencies' : """ filter out circular data dependencies [Gaudi::Algorithm] """,
    'InputCollection' : """ Input SimCalorimeterHit collection [BasicDigitizer] """,
    'OutputCollection' : """ Output SimCalorimeterHit collection [BasicDigitizer] """,
    'Threshold' : """ Minimum energy to keep a hit [BasicDigitizer] """,
    'DebugFrequency' : """ Print per-hit debug info every N events [BasicDigitizer] """,
  }
  __declaration_location__ = 'BasicDigitizer.cpp:97'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(BasicDigitizer, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'SND_reco'
  def getType( self ):
      return 'BasicDigitizer'
  pass # class BasicDigitizer

class ContributionInspector( ConfigurableAlgorithm ) :
  __slots__ = { 
    'ExtraInputs' : set(),
    'ExtraOutputs' : set(),
    'OutputLevel' : 0,
    'Enable' : True,
    'ErrorMax' : 1,
    'AuditInitialize' : False,
    'AuditReinitialize' : False,
    'AuditRestart' : False,
    'AuditExecute' : False,
    'AuditFinalize' : False,
    'AuditStart' : False,
    'AuditStop' : False,
    'Timeline' : True,
    'MonitorService' : 'MonitorSvc',
    'RegisterForContextService' : False,
    'Cardinality' : 0,
    'NeededResources' : [  ],
    'Asynchronous' : False,
    'FilterCircularDependencies' : True,
    'SiTargetCollection' : 'SiTargetHits',
    'SiPadCollection' : 'SiPadHits',
    'MaxHitsToPrint' : 20,
  }
  _propertyDocDct = { 
    'ExtraInputs' : """  [DataHandleHolderBase<PropertyHolder<CommonMessaging<implements<IAlgorithm,IDataHandleHolder,IProperty,IStateful> > > >] """,
    'ExtraOutputs' : """  [DataHandleHolderBase<PropertyHolder<CommonMessaging<implements<IAlgorithm,IDataHandleHolder,IProperty,IStateful> > > >] """,
    'OutputLevel' : """ output level [Gaudi::Algorithm] """,
    'Enable' : """ should the algorithm be executed or not [Gaudi::Algorithm] """,
    'ErrorMax' : """ [[deprecated]] max number of errors [Gaudi::Algorithm] """,
    'AuditInitialize' : """ trigger auditor on initialize() [Gaudi::Algorithm] """,
    'AuditReinitialize' : """ trigger auditor on reinitialize() [Gaudi::Algorithm] """,
    'AuditRestart' : """ trigger auditor on restart() [Gaudi::Algorithm] """,
    'AuditExecute' : """ trigger auditor on execute() [Gaudi::Algorithm] """,
    'AuditFinalize' : """ trigger auditor on finalize() [Gaudi::Algorithm] """,
    'AuditStart' : """ trigger auditor on start() [Gaudi::Algorithm] """,
    'AuditStop' : """ trigger auditor on stop() [Gaudi::Algorithm] """,
    'Timeline' : """ send events to TimelineSvc [Gaudi::Algorithm] """,
    'MonitorService' : """ name to use for Monitor Service [Gaudi::Algorithm] """,
    'RegisterForContextService' : """ flag to enforce the registration for Algorithm Context Service [Gaudi::Algorithm] """,
    'Cardinality' : """ how many clones to create - 0 means algo is reentrant [Gaudi::Algorithm] """,
    'NeededResources' : """ named resources needed during event looping [Gaudi::Algorithm] """,
    'Asynchronous' : """ whether algorithm is asynchronous and uses Boost Fiber to suspend while offloaded code is running. [Gaudi::Algorithm] """,
    'FilterCircularDependencies' : """ filter out circular data dependencies [Gaudi::Algorithm] """,
    'SiTargetCollection' : """ Input SiTarget SimCalorimeterHit collection [ContributionInspector] """,
    'SiPadCollection' : """ Input SiPad SimCalorimeterHit collection [ContributionInspector] """,
    'MaxHitsToPrint' : """ Max number of hits per collection for which contributions are printed at DEBUG level [ContributionInspector] """,
  }
  __declaration_location__ = 'ContributionInspector.cpp:118'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(ContributionInspector, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'SND_reco'
  def getType( self ):
      return 'ContributionInspector'
  pass # class ContributionInspector

class DelayTagger( ConfigurableAlgorithm ) :
  __slots__ = { 
    'ExtraInputs' : set(),
    'ExtraOutputs' : set(),
    'OutputLevel' : 0,
    'Enable' : True,
    'ErrorMax' : 1,
    'AuditInitialize' : False,
    'AuditReinitialize' : False,
    'AuditRestart' : False,
    'AuditExecute' : False,
    'AuditFinalize' : False,
    'AuditStart' : False,
    'AuditStop' : False,
    'Timeline' : True,
    'MonitorService' : 'MonitorSvc',
    'RegisterForContextService' : False,
    'Cardinality' : 0,
    'NeededResources' : [  ],
    'Asynchronous' : False,
    'FilterCircularDependencies' : True,
    'InputCollectionSiTarget' : 'SiTargetHits',
    'InputCollectionSiPad' : 'SiPadHits',
    'OutputCollectionSiTarget' : 'SiTargetHitsDelayed',
    'OutputCollectionSiPad' : 'SiPadHitsDelayed',
    'EventDelay' : 25.000000,
    'DebugFrequency' : 500,
  }
  _propertyDocDct = { 
    'ExtraInputs' : """  [DataHandleHolderBase<PropertyHolder<CommonMessaging<implements<IAlgorithm,IDataHandleHolder,IProperty,IStateful> > > >] """,
    'ExtraOutputs' : """  [DataHandleHolderBase<PropertyHolder<CommonMessaging<implements<IAlgorithm,IDataHandleHolder,IProperty,IStateful> > > >] """,
    'OutputLevel' : """ output level [Gaudi::Algorithm] """,
    'Enable' : """ should the algorithm be executed or not [Gaudi::Algorithm] """,
    'ErrorMax' : """ [[deprecated]] max number of errors [Gaudi::Algorithm] """,
    'AuditInitialize' : """ trigger auditor on initialize() [Gaudi::Algorithm] """,
    'AuditReinitialize' : """ trigger auditor on reinitialize() [Gaudi::Algorithm] """,
    'AuditRestart' : """ trigger auditor on restart() [Gaudi::Algorithm] """,
    'AuditExecute' : """ trigger auditor on execute() [Gaudi::Algorithm] """,
    'AuditFinalize' : """ trigger auditor on finalize() [Gaudi::Algorithm] """,
    'AuditStart' : """ trigger auditor on start() [Gaudi::Algorithm] """,
    'AuditStop' : """ trigger auditor on stop() [Gaudi::Algorithm] """,
    'Timeline' : """ send events to TimelineSvc [Gaudi::Algorithm] """,
    'MonitorService' : """ name to use for Monitor Service [Gaudi::Algorithm] """,
    'RegisterForContextService' : """ flag to enforce the registration for Algorithm Context Service [Gaudi::Algorithm] """,
    'Cardinality' : """ how many clones to create - 0 means algo is reentrant [Gaudi::Algorithm] """,
    'NeededResources' : """ named resources needed during event looping [Gaudi::Algorithm] """,
    'Asynchronous' : """ whether algorithm is asynchronous and uses Boost Fiber to suspend while offloaded code is running. [Gaudi::Algorithm] """,
    'FilterCircularDependencies' : """ filter out circular data dependencies [Gaudi::Algorithm] """,
    'InputCollectionSiTarget' : """ Input SiTarget SimCalorimeterHit collection [DelayTagger] """,
    'InputCollectionSiPad' : """ Input SiPad SimCalorimeterHit collection [DelayTagger] """,
    'OutputCollectionSiTarget' : """ Output SiTarget SimCalorimeterHit collection [DelayTagger] """,
    'OutputCollectionSiPad' : """ Output SiPad SimCalorimeterHit collection [DelayTagger] """,
    'EventDelay' : """ Time delay between consecutive input events (ns) [DelayTagger] """,
    'DebugFrequency' : """ Print debug/info lines every N events [DelayTagger] """,
  }
  __declaration_location__ = 'DelayTagger.cpp:135'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(DelayTagger, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'SND_reco'
  def getType( self ):
      return 'DelayTagger'
  pass # class DelayTagger

class EDM4HEP2RNTuple( ConfigurableAlgorithm ) :
  __slots__ = { 
    'ExtraInputs' : set(),
    'ExtraOutputs' : set(),
    'OutputLevel' : 0,
    'Enable' : True,
    'ErrorMax' : 1,
    'AuditInitialize' : False,
    'AuditReinitialize' : False,
    'AuditRestart' : False,
    'AuditExecute' : False,
    'AuditFinalize' : False,
    'AuditStart' : False,
    'AuditStop' : False,
    'Timeline' : True,
    'MonitorService' : 'MonitorSvc',
    'RegisterForContextService' : False,
    'Cardinality' : 0,
    'NeededResources' : [  ],
    'Asynchronous' : False,
    'FilterCircularDependencies' : True,
    'Collections' : [ 'SiTargetHitsWindowed' , 'SiPadHitsWindowed' ],
    'BitFields' : [ 'system:8,layer:8,slice:4,plane:1,strip:14' , 'system:8,layer:8,slice:4,x:9,y:9' ],
    'SourceIDParams' : [ 'SiTargetSourceIDs' , 'SiPadSourceIDs' ],
    'DetectorIDs' : [ 0 , 1 ],
    'NTupleNames' : [ 'SiTarget' , 'SiPad' ],
    'InputFile' : 'timewindows.root',
    'OutputFile' : 'ShipHits.root',
  }
  _propertyDocDct = { 
    'ExtraInputs' : """  [DataHandleHolderBase<PropertyHolder<CommonMessaging<implements<IAlgorithm,IDataHandleHolder,IProperty,IStateful> > > >] """,
    'ExtraOutputs' : """  [DataHandleHolderBase<PropertyHolder<CommonMessaging<implements<IAlgorithm,IDataHandleHolder,IProperty,IStateful> > > >] """,
    'OutputLevel' : """ output level [Gaudi::Algorithm] """,
    'Enable' : """ should the algorithm be executed or not [Gaudi::Algorithm] """,
    'ErrorMax' : """ [[deprecated]] max number of errors [Gaudi::Algorithm] """,
    'AuditInitialize' : """ trigger auditor on initialize() [Gaudi::Algorithm] """,
    'AuditReinitialize' : """ trigger auditor on reinitialize() [Gaudi::Algorithm] """,
    'AuditRestart' : """ trigger auditor on restart() [Gaudi::Algorithm] """,
    'AuditExecute' : """ trigger auditor on execute() [Gaudi::Algorithm] """,
    'AuditFinalize' : """ trigger auditor on finalize() [Gaudi::Algorithm] """,
    'AuditStart' : """ trigger auditor on start() [Gaudi::Algorithm] """,
    'AuditStop' : """ trigger auditor on stop() [Gaudi::Algorithm] """,
    'Timeline' : """ send events to TimelineSvc [Gaudi::Algorithm] """,
    'MonitorService' : """ name to use for Monitor Service [Gaudi::Algorithm] """,
    'RegisterForContextService' : """ flag to enforce the registration for Algorithm Context Service [Gaudi::Algorithm] """,
    'Cardinality' : """ how many clones to create - 0 means algo is reentrant [Gaudi::Algorithm] """,
    'NeededResources' : """ named resources needed during event looping [Gaudi::Algorithm] """,
    'Asynchronous' : """ whether algorithm is asynchronous and uses Boost Fiber to suspend while offloaded code is running. [Gaudi::Algorithm] """,
    'FilterCircularDependencies' : """ filter out circular data dependencies [Gaudi::Algorithm] """,
    'Collections' : """ Input SimCalorimeterHit collection names, one per detector [EDM4HEP2RNTuple] """,
    'BitFields' : """ DD4hep BitField encoding string for each collection [EDM4HEP2RNTuple] """,
    'SourceIDParams' : """ Frame parameter name carrying source_id vector for each collection [EDM4HEP2RNTuple] """,
    'DetectorIDs' : """ detector_id value (for logging) for each collection [EDM4HEP2RNTuple] """,
    'NTupleNames' : """ RNTuple name inside the output file for each collection [EDM4HEP2RNTuple] """,
    'InputFile' : """ Input edm4hep file (same file read by IOSvc) [EDM4HEP2RNTuple] """,
    'OutputFile' : """ Output ROOT file containing all per-detector RNTuples [EDM4HEP2RNTuple] """,
  }
  __declaration_location__ = 'EDM4HEP2RNTuple.cpp:360'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(EDM4HEP2RNTuple, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'SND_reco'
  def getType( self ):
      return 'EDM4HEP2RNTuple'
  pass # class EDM4HEP2RNTuple

class EventMerger( ConfigurableAlgorithm ) :
  __slots__ = { 
    'ExtraInputs' : set(),
    'ExtraOutputs' : set(),
    'OutputLevel' : 0,
    'Enable' : True,
    'ErrorMax' : 1,
    'AuditInitialize' : False,
    'AuditReinitialize' : False,
    'AuditRestart' : False,
    'AuditExecute' : False,
    'AuditFinalize' : False,
    'AuditStart' : False,
    'AuditStop' : False,
    'Timeline' : True,
    'MonitorService' : 'MonitorSvc',
    'RegisterForContextService' : False,
    'Cardinality' : 0,
    'NeededResources' : [  ],
    'Asynchronous' : False,
    'FilterCircularDependencies' : True,
    'InputCollectionSiTarget' : 'SiTargetHitsDelayed',
    'InputCollectionSiPad' : 'SiPadHitsDelayed',
    'OutputFile' : 'superevent.edm4hep.root',
    'OutputCollectionSiTarget' : 'SiTargetHitsMerged',
    'OutputCollectionSiPad' : 'SiPadHitsMerged',
  }
  _propertyDocDct = { 
    'ExtraInputs' : """  [DataHandleHolderBase<PropertyHolder<CommonMessaging<implements<IAlgorithm,IDataHandleHolder,IProperty,IStateful> > > >] """,
    'ExtraOutputs' : """  [DataHandleHolderBase<PropertyHolder<CommonMessaging<implements<IAlgorithm,IDataHandleHolder,IProperty,IStateful> > > >] """,
    'OutputLevel' : """ output level [Gaudi::Algorithm] """,
    'Enable' : """ should the algorithm be executed or not [Gaudi::Algorithm] """,
    'ErrorMax' : """ [[deprecated]] max number of errors [Gaudi::Algorithm] """,
    'AuditInitialize' : """ trigger auditor on initialize() [Gaudi::Algorithm] """,
    'AuditReinitialize' : """ trigger auditor on reinitialize() [Gaudi::Algorithm] """,
    'AuditRestart' : """ trigger auditor on restart() [Gaudi::Algorithm] """,
    'AuditExecute' : """ trigger auditor on execute() [Gaudi::Algorithm] """,
    'AuditFinalize' : """ trigger auditor on finalize() [Gaudi::Algorithm] """,
    'AuditStart' : """ trigger auditor on start() [Gaudi::Algorithm] """,
    'AuditStop' : """ trigger auditor on stop() [Gaudi::Algorithm] """,
    'Timeline' : """ send events to TimelineSvc [Gaudi::Algorithm] """,
    'MonitorService' : """ name to use for Monitor Service [Gaudi::Algorithm] """,
    'RegisterForContextService' : """ flag to enforce the registration for Algorithm Context Service [Gaudi::Algorithm] """,
    'Cardinality' : """ how many clones to create - 0 means algo is reentrant [Gaudi::Algorithm] """,
    'NeededResources' : """ named resources needed during event looping [Gaudi::Algorithm] """,
    'Asynchronous' : """ whether algorithm is asynchronous and uses Boost Fiber to suspend while offloaded code is running. [Gaudi::Algorithm] """,
    'FilterCircularDependencies' : """ filter out circular data dependencies [Gaudi::Algorithm] """,
    'InputCollectionSiTarget' : """ Input SiTarget collection [EventMerger] """,
    'InputCollectionSiPad' : """ Input SiPad collection [EventMerger] """,
    'OutputFile' : """ Output ROOT file for the super-event [EventMerger] """,
    'OutputCollectionSiTarget' : """ Output SiTarget collection name [EventMerger] """,
    'OutputCollectionSiPad' : """ Output SiPad collection name [EventMerger] """,
  }
  __declaration_location__ = 'EventMerger.cpp:149'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(EventMerger, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'SND_reco'
  def getType( self ):
      return 'EventMerger'
  pass # class EventMerger

class EventShuffler( ConfigurableAlgorithm ) :
  __slots__ = { 
    'ExtraInputs' : set(),
    'ExtraOutputs' : set(),
    'OutputLevel' : 0,
    'Enable' : True,
    'ErrorMax' : 1,
    'AuditInitialize' : False,
    'AuditReinitialize' : False,
    'AuditRestart' : False,
    'AuditExecute' : False,
    'AuditFinalize' : False,
    'AuditStart' : False,
    'AuditStop' : False,
    'Timeline' : True,
    'MonitorService' : 'MonitorSvc',
    'RegisterForContextService' : False,
    'Cardinality' : 0,
    'NeededResources' : [  ],
    'Asynchronous' : False,
    'FilterCircularDependencies' : True,
    'InputFiles' : [  ],
    'SourceIDs' : [  ],
    'Delays' : [  ],
    'CollectionsSiTarget' : [  ],
    'CollectionsSiPad' : [  ],
    'OutputFile' : 'shuffled.root',
    'OutputCollectionSiTarget' : 'SiTargetHitsMerged',
    'OutputCollectionSiPad' : 'SiPadHitsMerged',
  }
  _propertyDocDct = { 
    'ExtraInputs' : """  [DataHandleHolderBase<PropertyHolder<CommonMessaging<implements<IAlgorithm,IDataHandleHolder,IProperty,IStateful> > > >] """,
    'ExtraOutputs' : """  [DataHandleHolderBase<PropertyHolder<CommonMessaging<implements<IAlgorithm,IDataHandleHolder,IProperty,IStateful> > > >] """,
    'OutputLevel' : """ output level [Gaudi::Algorithm] """,
    'Enable' : """ should the algorithm be executed or not [Gaudi::Algorithm] """,
    'ErrorMax' : """ [[deprecated]] max number of errors [Gaudi::Algorithm] """,
    'AuditInitialize' : """ trigger auditor on initialize() [Gaudi::Algorithm] """,
    'AuditReinitialize' : """ trigger auditor on reinitialize() [Gaudi::Algorithm] """,
    'AuditRestart' : """ trigger auditor on restart() [Gaudi::Algorithm] """,
    'AuditExecute' : """ trigger auditor on execute() [Gaudi::Algorithm] """,
    'AuditFinalize' : """ trigger auditor on finalize() [Gaudi::Algorithm] """,
    'AuditStart' : """ trigger auditor on start() [Gaudi::Algorithm] """,
    'AuditStop' : """ trigger auditor on stop() [Gaudi::Algorithm] """,
    'Timeline' : """ send events to TimelineSvc [Gaudi::Algorithm] """,
    'MonitorService' : """ name to use for Monitor Service [Gaudi::Algorithm] """,
    'RegisterForContextService' : """ flag to enforce the registration for Algorithm Context Service [Gaudi::Algorithm] """,
    'Cardinality' : """ how many clones to create - 0 means algo is reentrant [Gaudi::Algorithm] """,
    'NeededResources' : """ named resources needed during event looping [Gaudi::Algorithm] """,
    'Asynchronous' : """ whether algorithm is asynchronous and uses Boost Fiber to suspend while offloaded code is running. [Gaudi::Algorithm] """,
    'FilterCircularDependencies' : """ filter out circular data dependencies [Gaudi::Algorithm] """,
    'InputFiles' : """ List of input edm4hep ROOT files [EventShuffler] """,
    'SourceIDs' : """ Source ID for each input file (same order as InputFiles) [EventShuffler] """,
    'Delays' : """ Inter-event delay in ns for each input file [EventShuffler] """,
    'CollectionsSiTarget' : """ SiTarget collection name for each input file [EventShuffler] """,
    'CollectionsSiPad' : """ SiPad collection name for each input file [EventShuffler] """,
    'OutputFile' : """ Output super-event ROOT file [EventShuffler] """,
    'OutputCollectionSiTarget' : """ Output SiTarget collection name [EventShuffler] """,
    'OutputCollectionSiPad' : """ Output SiPad collection name [EventShuffler] """,
  }
  __declaration_location__ = 'EventShuffler.cpp:288'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(EventShuffler, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'SND_reco'
  def getType( self ):
      return 'EventShuffler'
  pass # class EventShuffler

class EventWindowSplitter( ConfigurableAlgorithm ) :
  __slots__ = { 
    'ExtraInputs' : set(),
    'ExtraOutputs' : set(),
    'OutputLevel' : 0,
    'Enable' : True,
    'ErrorMax' : 1,
    'AuditInitialize' : False,
    'AuditReinitialize' : False,
    'AuditRestart' : False,
    'AuditExecute' : False,
    'AuditFinalize' : False,
    'AuditStart' : False,
    'AuditStop' : False,
    'Timeline' : True,
    'MonitorService' : 'MonitorSvc',
    'RegisterForContextService' : False,
    'Cardinality' : 0,
    'NeededResources' : [  ],
    'Asynchronous' : False,
    'FilterCircularDependencies' : True,
    'InputFile' : 'shuffled.root',
    'InputCollectionSiTarget' : 'SiTargetHitsMerged',
    'InputCollectionSiPad' : 'SiPadHitsMerged',
    'OutputFile' : 'timewindows.edm4hep.root',
    'OutputCollectionSiTarget' : 'SiTargetHitsWindowed',
    'OutputCollectionSiPad' : 'SiPadHitsWindowed',
    'WindowSize' : 25.000000,
  }
  _propertyDocDct = { 
    'ExtraInputs' : """  [DataHandleHolderBase<PropertyHolder<CommonMessaging<implements<IAlgorithm,IDataHandleHolder,IProperty,IStateful> > > >] """,
    'ExtraOutputs' : """  [DataHandleHolderBase<PropertyHolder<CommonMessaging<implements<IAlgorithm,IDataHandleHolder,IProperty,IStateful> > > >] """,
    'OutputLevel' : """ output level [Gaudi::Algorithm] """,
    'Enable' : """ should the algorithm be executed or not [Gaudi::Algorithm] """,
    'ErrorMax' : """ [[deprecated]] max number of errors [Gaudi::Algorithm] """,
    'AuditInitialize' : """ trigger auditor on initialize() [Gaudi::Algorithm] """,
    'AuditReinitialize' : """ trigger auditor on reinitialize() [Gaudi::Algorithm] """,
    'AuditRestart' : """ trigger auditor on restart() [Gaudi::Algorithm] """,
    'AuditExecute' : """ trigger auditor on execute() [Gaudi::Algorithm] """,
    'AuditFinalize' : """ trigger auditor on finalize() [Gaudi::Algorithm] """,
    'AuditStart' : """ trigger auditor on start() [Gaudi::Algorithm] """,
    'AuditStop' : """ trigger auditor on stop() [Gaudi::Algorithm] """,
    'Timeline' : """ send events to TimelineSvc [Gaudi::Algorithm] """,
    'MonitorService' : """ name to use for Monitor Service [Gaudi::Algorithm] """,
    'RegisterForContextService' : """ flag to enforce the registration for Algorithm Context Service [Gaudi::Algorithm] """,
    'Cardinality' : """ how many clones to create - 0 means algo is reentrant [Gaudi::Algorithm] """,
    'NeededResources' : """ named resources needed during event looping [Gaudi::Algorithm] """,
    'Asynchronous' : """ whether algorithm is asynchronous and uses Boost Fiber to suspend while offloaded code is running. [Gaudi::Algorithm] """,
    'FilterCircularDependencies' : """ filter out circular data dependencies [Gaudi::Algorithm] """,
    'InputFile' : """ Input super-event file [EventWindowSplitter] """,
    'InputCollectionSiTarget' : """ Input SiTarget collection [EventWindowSplitter] """,
    'InputCollectionSiPad' : """ Input SiPad collection [EventWindowSplitter] """,
    'OutputFile' : """ Output ROOT file for time-windowed events [EventWindowSplitter] """,
    'OutputCollectionSiTarget' : """ Output SiTarget collection name [EventWindowSplitter] """,
    'OutputCollectionSiPad' : """ Output SiPad collection name [EventWindowSplitter] """,
    'WindowSize' : """ Time window size (ns) [EventWindowSplitter] """,
  }
  __declaration_location__ = 'EventWindowSplitter.cpp:290'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(EventWindowSplitter, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'SND_reco'
  def getType( self ):
      return 'EventWindowSplitter'
  pass # class EventWindowSplitter

class GeV2MIPConversion( ConfigurableAlgorithm ) :
  __slots__ = { 
    'ExtraInputs' : set(),
    'ExtraOutputs' : set(),
    'OutputLevel' : 0,
    'Enable' : True,
    'ErrorMax' : 1,
    'AuditInitialize' : False,
    'AuditReinitialize' : False,
    'AuditRestart' : False,
    'AuditExecute' : False,
    'AuditFinalize' : False,
    'AuditStart' : False,
    'AuditStop' : False,
    'Timeline' : True,
    'MonitorService' : 'MonitorSvc',
    'RegisterForContextService' : False,
    'Cardinality' : 0,
    'NeededResources' : [  ],
    'Asynchronous' : False,
    'FilterCircularDependencies' : True,
    'InputCollection' : 'SiTargetHits',
    'OutputCollection' : 'SiTargetHitsMIP',
    'MIPValue' : 9.0000000e-06,
    'DebugFrequency' : 500,
  }
  _propertyDocDct = { 
    'ExtraInputs' : """  [DataHandleHolderBase<PropertyHolder<CommonMessaging<implements<IAlgorithm,IDataHandleHolder,IProperty,IStateful> > > >] """,
    'ExtraOutputs' : """  [DataHandleHolderBase<PropertyHolder<CommonMessaging<implements<IAlgorithm,IDataHandleHolder,IProperty,IStateful> > > >] """,
    'OutputLevel' : """ output level [Gaudi::Algorithm] """,
    'Enable' : """ should the algorithm be executed or not [Gaudi::Algorithm] """,
    'ErrorMax' : """ [[deprecated]] max number of errors [Gaudi::Algorithm] """,
    'AuditInitialize' : """ trigger auditor on initialize() [Gaudi::Algorithm] """,
    'AuditReinitialize' : """ trigger auditor on reinitialize() [Gaudi::Algorithm] """,
    'AuditRestart' : """ trigger auditor on restart() [Gaudi::Algorithm] """,
    'AuditExecute' : """ trigger auditor on execute() [Gaudi::Algorithm] """,
    'AuditFinalize' : """ trigger auditor on finalize() [Gaudi::Algorithm] """,
    'AuditStart' : """ trigger auditor on start() [Gaudi::Algorithm] """,
    'AuditStop' : """ trigger auditor on stop() [Gaudi::Algorithm] """,
    'Timeline' : """ send events to TimelineSvc [Gaudi::Algorithm] """,
    'MonitorService' : """ name to use for Monitor Service [Gaudi::Algorithm] """,
    'RegisterForContextService' : """ flag to enforce the registration for Algorithm Context Service [Gaudi::Algorithm] """,
    'Cardinality' : """ how many clones to create - 0 means algo is reentrant [Gaudi::Algorithm] """,
    'NeededResources' : """ named resources needed during event looping [Gaudi::Algorithm] """,
    'Asynchronous' : """ whether algorithm is asynchronous and uses Boost Fiber to suspend while offloaded code is running. [Gaudi::Algorithm] """,
    'FilterCircularDependencies' : """ filter out circular data dependencies [Gaudi::Algorithm] """,
    'InputCollection' : """ Input SimCalorimeterHit collection [GeV2MIPConversion] """,
    'OutputCollection' : """ Output SimCalorimeterHit collection [GeV2MIPConversion] """,
    'MIPValue' : """ MIP value for conversion (in GeV) [GeV2MIPConversion] """,
    'DebugFrequency' : """ Print per-hit debug info every N events [GeV2MIPConversion] """,
  }
  __declaration_location__ = 'GeV2MIPConversion.cpp:109'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(GeV2MIPConversion, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'SND_reco'
  def getType( self ):
      return 'GeV2MIPConversion'
  pass # class GeV2MIPConversion
