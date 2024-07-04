package archie

import archie.domain.pipelines.*
import archie.domain.pipelines.PipelineAssetsCollector
import archie.domain.pipelines.PipelineNamespace

import bpipe.PipelineEvent

import groovy.json.*

import java.nio.file.Path
import java.nio.file.Paths
import java.time.format.DateTimeFormatter
import java.time.OffsetDateTime
import java.time.ZonedDateTime
import java.time.ZoneId
import java.util.regex.Pattern

import static archie.apiclient.dto.AnalysisStatus.FAILED
import static archie.apiclient.dto.AnalysisStatus.IN_PROGRESS
import static archie.apiclient.dto.AnalysisStatus.SUCCESS
import static archie.domain.pipelines.PipelineStageName.PIPELINE_FINISHED
import static archie.utils.AnalysisUtils.createAnalysisUpdateMsgData
import static archie.Utils.sendMessage


class ArchieEventHandlers {

    static def onPipelineStarted = { String analysis_id, String project, PipelineEvent type, String desc, Map<String, Object> details ->
        def pipeline_id = bpipe.Config.config.pid
        println "Pipeline $type event triggered, pipeline_id=$pipeline_id, desc=$desc"

        def msgContentData = createAnalysisUpdateMsgData(analysis_id, project, pipeline_id, IN_PROGRESS.name(), Collections.emptyList(), null, IN_PROGRESS.name())

        sendMessage(msgContentData)
    }

    static def onPipelineFinished = { String analysis_id, String project, def meta, PipelineAssetsCollector assetsCollector, PipelineEvent type, String desc, Map<String, Object> details ->
        def pipeline_id = bpipe.Config.config.pid
        println "Pipeline $type event triggered, pipeline_id=$pipeline_id, desc=$desc"

        if(details.result) {
            Closure<PipelineStageAssets> stageAssetCollector = assetsCollector?.findStageAssetCollector(PIPELINE_FINISHED.id)
            if (stageAssetCollector) {
                Path analysisDir = Paths.get(new File('.').canonicalPath)
                PipelineStageAssets stageAssets = stageAssetCollector(meta, analysisDir)
                def msgContentData = createAnalysisUpdateMsgData(
                    analysis_id,
                    project,
                    pipeline_id,
                    SUCCESS.name(),
                    stageAssets.getAssets(),
                    PIPELINE_FINISHED.id,
                    desc,
                    stageAssets.getMetadata()
                )
                sendMessage(msgContentData)
            }
            else {
                def msgContentData = createAnalysisUpdateMsgData(analysis_id, project, pipeline_id, SUCCESS.name())
                sendMessage(msgContentData)
            }
        }
        else {
            println "Pipeline failed"
            def metadata = ['error': desc]
            def msgContentData = createAnalysisUpdateMsgData(analysis_id, project, pipeline_id, FAILED.name())
            msgContentData.pipeline_status = FAILED.name()
            msgContentData.metadata = metadata
            sendMessage(msgContentData)
        }
    }

    static void init_hook(Object meta, String analysis_id, String project, String pipeline_script) {
        println "init_hook: Samples are $meta"
        assert "$analysis_id", "analysis_id parameter must be defined"
        assert "$project", "project parameter must be defined"
        assert "$pipeline_script", "pipeline_script parameter must be defined"

        def sampleIdentifiers = meta*.key
        assert sampleIdentifiers && !sampleIdentifiers.empty, "sampleIdentifiers is invalid, make sure samples_parser is correctly defined and returns a map keyed by sample identifiers."

        String scriptName = Paths.get("$pipeline_script").fileName
        def assetsCollector = PipelineAssetsCollector.find(scriptName, PipelineNamespace.WARPY)

        println "init_hook: Adding Archie pipeline handlers for analysis_id=$analysis_id, project=$project, pipeline_script=$pipeline_script"
        bpipe.EventManager.instance.addListener(bpipe.PipelineEvent.STARTED, ArchieEventHandlers.onPipelineStarted.curry("$analysis_id", "$project"))
        bpipe.EventManager.instance.addListener(bpipe.PipelineEvent.FINISHED, ArchieEventHandlers.onPipelineFinished.curry("$analysis_id", "$project", meta, assetsCollector))
    }
}
