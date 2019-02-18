c = get_config()

c.Exporter.template_file = '/templates/latex.tplx'
c.TagRemovePreprocessor.remove_input_tags.add("rm_in")
c.TagRemovePreprocessor.remove_all_outputs_tags.add("rm_out")

# the following does the equivalent of --no-prompt, see here: https://github.com/jupyter/nbconvert/blob/master/nbconvert/nbconvertapp.py#L109-L114
# exporter_settings = {
#     'exclude_input_prompt' : True,
#     'exclude_output_prompt' : True,
# }
# c.TemplateExporter.update(exporter_settings)
