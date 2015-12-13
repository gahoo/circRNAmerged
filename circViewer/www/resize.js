$(function () {
  $(".resizable").resizable({ stop: function () {
    $(".shiny-plot-output").height($(this).height()*0.8);
    $(this).height('auto');
    $(this).css({"position":"fixed"});
    if( $(this).height() > $( window ).height() ){
      $(this).height($( window ).height());
      $(this).css({"overflow-y":"auto"});
    }else{
      $(this).css({"overflow-y":""});
      $(this).height('auto');
      $(this).width('auto');
    }
  } 
  });
});